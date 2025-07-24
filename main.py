from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import pandas as pd
import re
import os
from itertools import product
from functools import lru_cache
from itertools import product

fasta_path = "./data/avengers-3.fa"
primer_type = "RACE"  # Options: "qPCR" or "RACE"
mode_text = "coverage"  # Options: "specificity" or "coverage"

params_qpcr = {
    "fw_start": 50, 
    "fw_end": 100,
    "rev_start": 150, 
    "rev_end": 200,
    "k_range": range(18, 26),
    "gc_min": 50, 
    "gc_max": 70,
    "prod_min": 80, 
    "prod_max": 150,
    "tm_min": 57, 
    "tm_max": 62,
}

params_race = {
    "window_size": 200,
    "k_range": range(23, 29),
    "gc_min": 50, 
    "gc_max": 70,
    "tm_min": 57,
    "tm_max": 62,
}

def calcular_gc(seq):
    """
    Calculates the GC content of a sequence.
    """
    if not seq:
        return 0.0
    gc_count = seq.count("G") + seq.count("C")
    return round(gc_count / len(seq) * 100, 2)

def extract_kmers(seq, k_range):
    """
    Extracts all k-mers from a sequence for a given range of k values.
    """
    if not seq:
        return
    for k in k_range:
        if k > 0 and k <= len(seq):
            for i in range(len(seq) - k + 1):
                yield seq[i:i+k]

COMPLEMENT_MAP = {"A": "T", "T": "A", "G": "C", "C": "G"}

def is_complement(base1, base2):
    """
    Checks if two bases are complementary.
    """
    return COMPLEMENT_MAP.get(base1) == base2

def reverse_complement(seq):
    """
    Computes the reverse complement of a sequence.
    """
    return str(Seq(seq).reverse_complement())

def longest_complementary(s, t):
    """
    Finds the length of the longest substring in 's' that is perfectly
    complementary to a substring in 't'.
    """
    m, n = len(s), len(t)
    if m == 0 or n == 0:
        return 0

    dp = [0] * (n + 1)
    max_len = 0

    for i in range(1, m + 1):
        prev_row_val = 0
        for j in range(1, n + 1):
            current_dp_j = dp[j]
            if is_complement(s[i-1], t[j-1]):
                dp[j] = prev_row_val + 1
                if dp[j] > max_len:
                    max_len = dp[j]
            else:
                dp[j] = 0
            prev_row_val = current_dp_j
    return max_len

def calcular_tm(seq, Na=50e-3):
    """
    Calculates the melting temperature (Tm).
    """
    try:
        return round(mt.Tm_NN(str(seq), Na=Na), 1)
    except Exception:
        return 0.0

def extract_isoform(nombre):
    """
    Extracts an isoform number from a name using a regular expression.
    """
    m = re.search(r"-([0-9]+)", nombre)
    return m.group(1) if m else ""

def verificar_primer(primer, seqs):
    """
    Finds all sequence records that contain a given primer.
    """
    hits = [rec.id for rec in seqs if primer in str(rec.seq)]
    return {"total": len(hits), "isoformas": ",".join(hits)}

def get_candidates(region, k_range, gc_min, gc_max):
    """Finds unique k-mers in a region that fit GC content criteria."""
    candidates = set()
    for kmer in extract_kmers(region, k_range):
        if gc_min <= calcular_gc(kmer) <= gc_max:
            candidates.add(kmer)
    return list(candidates)

print(f"📂 Reading FASTA {fasta_path}")
seq_set = list(SeqIO.parse(fasta_path, "fasta"))
print(f"✅ {len(seq_set)} sequences loaded.")


@lru_cache(maxsize=1000)
def check_cached_primer(primer, seq_tuple):
    """Cached version of verificar_primer to avoid redundant computations"""
    seqs = list(seq_tuple)
    return verificar_primer(primer, seqs)

def prepare_sequences_for_cache(seqs):
    """Convert sequences to a hashable format for caching"""
    return tuple((rec.id, str(rec.seq)) for rec in seqs)

def generar_qPCR(seqs, params, mode):
    """Optimized qPCR primer generation with reduced loops and caching"""
    seq_tuple = prepare_sequences_for_cache(seqs)
    resb = {}
    reso = {}
    
    valid_seqs = [(rec.id, str(rec.seq)) for rec in seqs 
                  if len(str(rec.seq)) >= params["rev_end"]]
    
    if not valid_seqs:
        return {"best": resb, "other": reso}
    
    for name, seq in valid_seqs:
        fwR = seq[params["fw_start"] - 1:params["fw_end"]]
        rvR = seq[params["rev_start"] - 1:params["rev_end"]]
        
        fw_cands = get_candidates(fwR, params["k_range"], params["gc_min"], params["gc_max"])
        rv_cands = get_candidates(rvR, params["k_range"], params["gc_min"], params["gc_max"])
        
        if not fw_cands or not rv_cands:
            resb[name] = None
            reso[name] = []
            continue
        
        # Pre-verify all primers to avoid redundant calls
        fw_checks = {fw: check_cached_primer(fw, seq_tuple) for fw in fw_cands}
        rv_checks = {rv: check_cached_primer(rv, seq_tuple) for rv in rv_cands}
        
        # Pre-filter valid primers for specificity mode
        if mode == "specificity":
            iso = extract_isoform(name)
            valid_fw = [fw for fw, chk in fw_checks.items() 
                       if chk["total"] == 1 and iso in chk["isoformas"]]
            valid_rv = [rv for rv, chk in rv_checks.items() 
                       if chk["total"] == 1 and iso in chk["isoformas"]]
        else:
            valid_fw = fw_cands
            valid_rv = rv_cands
        
        best = None
        others = []
        
        for fw, rv in product(valid_fw, valid_rv):
            size = params["rev_end"] - params["fw_start"] + len(rv)
            if not (params["prod_min"] <= size <= params["prod_max"]):
                continue
                
            primer = {"forward": fw, "reverse": rv, "product_size": size}
            
            if mode == "specificity":
                # Already filtered for valid primers above
                best = primer
                break
            else:
                if len(others) < 10:
                    others.append(primer)
                if not best:  # Take first as fallback
                    best = primer
        
        resb[name] = best
        reso[name] = others
    
    return {"best": resb, "other": reso}

def generar_RACE(seqs, params, mode):
    """Optimized RACE primer generation"""
    seq_tuple = prepare_sequences_for_cache(seqs)
    resb = {}
    reso = {}
    
    # Pre-filter sequences by length
    valid_seqs = [(rec.id, str(rec.seq)) for rec in seqs 
                  if len(str(rec.seq)) >= 2 * params["window_size"]]
    
    if not valid_seqs:
        return {"best": resb, "other": reso}
    
    for name, seq in valid_seqs:
        # Extract regions once
        fwR = seq[:params["window_size"]]
        rvR = seq[-params["window_size"]:]
        
        # Get candidates once
        fw_cands = get_candidates(fwR, params["k_range"], params["gc_min"], params["gc_max"])
        rv_cands = get_candidates(rvR, params["k_range"], params["gc_min"], params["gc_max"])
        
        if not fw_cands or not rv_cands:
            resb[name] = None
            reso[name] = []
            continue
        
        # Pre-verify all primers
        fw_checks = {fw: check_cached_primer(fw, seq_tuple) for fw in fw_cands}
        rv_checks = {rv: check_cached_primer(rv, seq_tuple) for rv in rv_cands}
        
        # Pre-filter for specificity mode
        if mode == "specificity":
            iso = extract_isoform(name)
            valid_fw = [fw for fw, chk in fw_checks.items() 
                       if chk["total"] == 1 and iso in chk["isoformas"]]
            valid_rv = [rv for rv, chk in rv_checks.items() 
                       if chk["total"] == 1 and iso in chk["isoformas"]]
        else:
            valid_fw = fw_cands
            valid_rv = rv_cands
        
        best = None
        others = []
        
        # Use product for cleaner iteration
        for fw, rv in product(valid_fw, valid_rv):
            primer = {"forward": fw, "reverse": rv, "product_size": None}
            
            if mode == "specificity":
                # Already filtered for valid primers
                best = primer
                break
            else:
                if len(others) < 10:
                    others.append(primer)
                if not best:
                    best = primer
        
        resb[name] = best
        reso[name] = others
    
    return {"best": resb, "other": reso}

def generar_coverage(seqs, params_qpcr, params_race, primer_type):
    """Optimized coverage generation using set operations"""
    iso_names = set(rec.id for rec in seqs)  # Use set for faster operations
    
    # Determine parameters and region extraction function
    if primer_type == "qPCR":
        params = params_qpcr
        get_regions = lambda s: (
            s[params["fw_start"] - 1:params["fw_end"]],
            s[params["rev_start"] - 1:params["rev_end"]]
        )
    else:
        params = params_race
        get_regions = lambda s: (
            s[:params["window_size"]],
            s[-params["window_size"]:]
        )
    
    # Pre-compute all k-mers and their hits in one pass
    fw_kmer_hits = {}
    rv_kmer_hits = {}
    
    for rec in seqs:
        seq_str = str(rec.seq)
        fw_region, rv_region = get_regions(seq_str)
        
        # Extract k-mers for forward region
        fw_kmers = extract_kmers(fw_region, params["k_range"])
        for kmer in fw_kmers:
            if params["gc_min"] <= calcular_gc(kmer) <= params["gc_max"]:
                if kmer not in fw_kmer_hits:
                    fw_kmer_hits[kmer] = set()
                fw_kmer_hits[kmer].add(rec.id)
        
        # Extract k-mers for reverse region
        rv_kmers = extract_kmers(rv_region, params["k_range"])
        for kmer in rv_kmers:
            if params["gc_min"] <= calcular_gc(kmer) <= params["gc_max"]:
                if kmer not in rv_kmer_hits:
                    rv_kmer_hits[kmer] = set()
                rv_kmer_hits[kmer].add(rec.id)
    
    # Generate primer combinations and find coverage
    selected = []
    pending = iso_names.copy()
    
    # Sort combinations by coverage potential (descending)
    combinations = []
    for fw_kmer, fw_hits in fw_kmer_hits.items():
        for rv_kmer, rv_hits in rv_kmer_hits.items():
            common_hits = fw_hits & rv_hits & pending
            if common_hits:  # Only consider combinations that cover pending isoforms
                combinations.append((len(common_hits), fw_kmer, rv_kmer, common_hits))
    
    # Sort by coverage size (descending) to prioritize high-coverage primers
    combinations.sort(reverse=True)
    
    for coverage_size, fw_kmer, rv_kmer, common_hits in combinations:
        if not pending:  # All isoforms covered
            break
        
        # Only add if this combination covers some pending isoforms
        actual_coverage = common_hits & pending
        if actual_coverage:
            selected.append({
                "forward": fw_kmer, 
                "reverse": rv_kmer, 
                "isoforms": list(actual_coverage)
            })
            pending -= actual_coverage
    
    return {"selected": selected, "pending": list(pending)}

if mode_text == "specificity":
    if primer_type == "qPCR":
        resultados = generar_qPCR(seq_set, params_qpcr, mode_text)
    elif primer_type == "RACE":
        resultados = generar_RACE(seq_set, params_race, mode_text)
    else:
        raise BaseException("⚠️ Unsupported primer type")
elif mode_text == "coverage":
    resultados = generar_coverage(seq_set, params_qpcr, params_race, primer_type)
else:
    raise BaseException("⚠️ Unsupported mode")

# Create output directory if it doesn't exist
os.makedirs("./output", exist_ok=True)

if mode_text == "specificity":
    rows = []
    for nm, p in resultados["best"].items():
        if not p:
            continue
        f, r = p["forward"], p["reverse"]
        row = {
            "Isoforma": nm,
            "Primer_Forward": f,
            "Primer_Reverse": r,
            "GC_Forward": calcular_gc(f),
            "GC_Reverse": calcular_gc(r),
            "Tm_Forward": calcular_tm(f),
            "Tm_Reverse": calcular_tm(r),
            "Hairpin_Forward": longest_complementary(f, reverse_complement(f)),
            "Hairpin_Reverse": longest_complementary(r, reverse_complement(r)),
            "SelfDimer_Forward": longest_complementary(f, f),
            "SelfDimer_Reverse": longest_complementary(r, r),
            "CrossDimer": longest_complementary(f, reverse_complement(r))
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    output_path = "./output/results_specificity.csv"

elif mode_text == "coverage":
    rows = []
    for item in resultados["selected"]:
        row = {
            "Primer_Forward": item["forward"],
            "Primer_Reverse": item["reverse"],
            "Isoformas": ";".join(item["isoforms"])
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    output_path = "./output/results_coverage.csv"

else:
    raise BaseException("⚠️ Unsupported mode")

df.to_csv(output_path, index=False)
print(f"✅ Results saved to {output_path}")