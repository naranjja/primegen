# --- Standard Library Imports ---
import re
import os
from itertools import product
from functools import lru_cache

# --- Third-party Imports ---
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

# ==============================================================================
# --- CONFIGURATION PARAMETERS ---
# ==============================================================================

# --- qPCR Specific Parameters ---
# Defines the regions to search for forward and reverse primers and product size constraints.
PARAMS_QPCR = {
    "fw_start": 50,             # Start position for forward primer search region.
    "fw_end": 100,              # End position for forward primer search region.
    "rev_start": 150,           # Start position for reverse primer search region.
    "rev_end": 200,             # End position for reverse primer search region.
    "k_range": range(18, 26),   # Range of possible k-mer (primer) lengths.
    "gc_min": 50,               # Minimum GC content percentage.
    "gc_max": 70,               # Maximum GC content percentage.
    "prod_min": 80,             # Minimum PCR product size.
    "prod_max": 150,            # Maximum PCR product size.
    "tm_min": 57,               # Minimum melting temperature (Tm).
    "tm_max": 62,               # Maximum melting temperature (Tm).
}

# --- RACE Specific Parameters ---
# Defines the search window at the ends of the sequence for RACE primers.
PARAMS_RACE = {
    "window_size": 200,         # Size of the search window from the 5' and 3' ends.
    "k_range": range(23, 29),   # Range of possible k-mer (primer) lengths.
    "gc_min": 50,               # Minimum GC content percentage.
    "gc_max": 70,               # Maximum GC content percentage.
    "tm_min": 57,               # Minimum melting temperature (Tm).
    "tm_max": 62,               # Maximum melting temperature (Tm).
}

# ==============================================================================
# --- CORE BIOINFORMATICS UTILITY FUNCTIONS ---
# ==============================================================================

def calculate_gc_content(sequence):
    """
    Calculates the GC content of a DNA sequence.

    Args:
        sequence (str): The DNA sequence string.

    Returns:
        float: The percentage of G and C bases, rounded to two decimal places.
               Returns 0.0 if the sequence is empty.
    """
    if not sequence:
        return 0.0
    gc_count = sequence.count("G") + sequence.count("C")
    return round(gc_count / len(sequence) * 100, 2)

def extract_kmers(sequence, k_range):
    """
    Yields all k-mers of specified lengths from a sequence.

    This is a generator function that produces k-mers on-the-fly to save memory,
    rather than storing them all in a list.

    Args:
        sequence (str): The DNA sequence.
        k_range (range): A range of k-mer lengths to extract (e.g., range(18, 26)).

    Yields:
        str: The next k-mer found in the sequence.
    """
    if not sequence:
        return
    for k in k_range:
        if 0 < k <= len(sequence):
            for i in range(len(sequence) - k + 1):
                yield sequence[i:i+k]

# Pre-computed map for fast base complementarity checks.
COMPLEMENT_MAP = {"A": "T", "T": "A", "G": "C", "C": "G"}

def are_bases_complementary(base1, base2):
    """Checks if two DNA bases are complementary using a fast lookup."""
    return COMPLEMENT_MAP.get(base1) == base2

def reverse_complement(sequence):
    """Computes the reverse complement of a DNA sequence using Biopython."""
    return str(Seq(sequence).reverse_complement())

def find_longest_complementary_run(seq1, seq2):
    """
    Finds the length of the longest perfectly complementary substring between two sequences.
    This is used to predict secondary structures like hairpins and dimers.

    This function uses a memory-optimized dynamic programming approach (O(n) space).

    Args:
        seq1 (str): The first DNA sequence.
        seq2 (str): The second DNA sequence.

    Returns:
        int: The length of the longest complementary run.
    """
    m, n = len(seq1), len(seq2)
    if m == 0 or n == 0:
        return 0

    dp_row = [0] * (n + 1)
    max_len = 0

    for i in range(1, m + 1):
        prev_diagonal_val = 0
        for j in range(1, n + 1):
            current_dp_j = dp_row[j]
            if are_bases_complementary(seq1[i-1], seq2[j-1]):
                dp_row[j] = prev_diagonal_val + 1
                if dp_row[j] > max_len:
                    max_len = dp_row[j]
            else:
                dp_row[j] = 0
            prev_diagonal_val = current_dp_j
    return max_len

def calculate_melting_temp(sequence, na_conc=50e-3):
    """
    Calculates the melting temperature (Tm) using Biopython's nearest-neighbor model.
    """
    try:
        return round(mt.Tm_NN(str(sequence), Na=na_conc), 1)
    except Exception:
        return 0.0

def extract_isoform_id(name):
    """
    Extracts an isoform number (e.g., '3') from a sequence name.
    Tries to match patterns like 'gene-3', 'transcript_2', etc.
    """
    match = re.search(r"[-_]?(\d+)$", name)
    if not match:
        print(f"⚠️ Could not extract isoform ID from: {name}")
    return match.group(1) if match else ""

def check_primer_hits(primer, sequences):
    """
    Finds all sequence records that contain a given primer sequence
    or its reverse complement.
    """
    primer_rc = reverse_complement(primer)
    hits = [
        rec.id
        for rec in sequences
        if primer in str(rec.seq) or primer_rc in str(rec.seq)
    ]
    return {"total": len(hits), "isoforms": ",".join(hits)}

def find_candidate_primers(region, k_range, gc_min, gc_max):
    """
    Finds unique k-mers within a sequence region that meet GC content criteria.

    Args:
        region (str): The DNA sequence region to search.
        k_range (range): The range of k-mer lengths.
        gc_min (float): The minimum allowed GC content.
        gc_max (float): The maximum allowed GC content.

    Returns:
        list: A list of unique candidate primer sequences.
    """
    candidates = set()
    for kmer in extract_kmers(region, k_range):
        if gc_min <= calculate_gc_content(kmer) <= gc_max:
            candidates.add(kmer)
    return list(candidates)

# ==============================================================================
# --- PRIMER DESIGN CORE LOGIC ---
# ==============================================================================

@lru_cache(maxsize=10000)
def check_cached_primer_hits(primer, sequences_tuple):
    """
    A cached wrapper for `check_primer_hits` to avoid re-computing hits
    for the same primer. The `lru_cache` decorator memoizes the results.
    """
    # The cache requires hashable arguments, so we convert the tuple back to a list of records.
    sequences = [SeqIO.SeqRecord(Seq(seq), id=name) for name, seq in sequences_tuple]
    return check_primer_hits(primer, sequences)

def prepare_sequences_for_caching(sequences):
    """Converts a list of SeqRecord objects to a hashable tuple for caching."""
    return tuple((rec.id, str(rec.seq)) for rec in sequences)

def design_qpcr_primers(sequences, params, mode):
    """
    Main logic for designing qPCR primers for either specificity or coverage.
    """
    sequences_tuple = prepare_sequences_for_caching(sequences)
    best_primers = {}
    other_options = {}

    valid_sequences = [(rec.id, str(rec.seq)) for rec in sequences
                       if len(str(rec.seq)) >= params["rev_end"]]

    if not valid_sequences:
        return {"best": best_primers, "other": other_options}

    for name, seq in valid_sequences:
        fw_region = seq[params["fw_start"] - 1:params["fw_end"]]
        rev_region = seq[params["rev_start"] - 1:params["rev_end"]]

        fw_candidates = find_candidate_primers(fw_region, params["k_range"], params["gc_min"], params["gc_max"])
        # Find k-mers on the template strand, then reverse-complement them to get actual primers.
        rev_template_candidates = find_candidate_primers(rev_region, params["k_range"], params["gc_min"], params["gc_max"])
        # Create a map of {actual_reverse_primer: template_sequence} for calculating product size later.
        rev_candidates_map = {reverse_complement(p): p for p in rev_template_candidates}

        if not fw_candidates or not rev_candidates_map:
            best_primers[name] = None
            other_options[name] = []
            continue

        fw_checks = {fw: check_cached_primer_hits(fw, sequences_tuple) for fw in fw_candidates}
        # Check specificity using the actual reverse primer sequences.
        rev_checks = {rv: check_cached_primer_hits(rv, sequences_tuple) for rv in rev_candidates_map.keys()}

        valid_fw = fw_candidates
        valid_rev = list(rev_candidates_map.keys())

        if mode == "specificity":
            isoform_id = extract_isoform_id(name)
            valid_fw = [
                fw for fw, check in fw_checks.items()
                if isoform_id in check["isoforms"] and check["total"] <= 2
            ]
            valid_rev = [
                rv for rv, check in rev_checks.items()
                if isoform_id in check["isoforms"] and check["total"] <= 2
            ]
        
        print(f"🔍 [{name}] Found {len(valid_fw)} FW primers, {len(valid_rev)} REV primers after specificity filtering.")

        best_pair = None
        other_pairs = []

        for fw, rv in product(valid_fw, valid_rev):
            # --- CORRECTED PRODUCT SIZE CALCULATION ---
            fw_pos = seq.find(fw, params["fw_start"] - 1)
            # Retrieve the original template sequence for the reverse primer to find its position.
            rv_template = rev_candidates_map[rv]
            rv_template_pos = seq.find(rv_template, params["rev_start"] - 1)

            if fw_pos == -1 or rv_template_pos == -1:
                continue

            # Correctly calculate product size.
            product_size = (rv_template_pos + len(rv_template)) - fw_pos
            
            if not (params["prod_min"] <= product_size <= params["prod_max"]):
                continue

            primer_pair = {"forward": fw, "reverse": rv, "product_size": product_size}

            if mode == "specificity":
                best_pair = primer_pair
                break  # Found the first specific pair, stop searching for this isoform.
            else: # coverage mode
                if len(other_pairs) < 10:
                    other_pairs.append(primer_pair)
                if not best_pair:
                    best_pair = primer_pair

        # --- REMOVED LOGIC FLAW ---
        # The original code had a 'break' here that prematurely exited the main loop.

        best_primers[name] = best_pair
        other_options[name] = other_pairs

    return {"best": best_primers, "other": other_options}

def design_race_primers(sequences, params, mode):
    """
    Main logic for designing RACE primers for either specificity or coverage.
    """
    sequences_tuple = prepare_sequences_for_caching(sequences)
    best_primers = {}
    other_options = {}

    valid_sequences = [(rec.id, str(rec.seq)) for rec in sequences
                       if len(str(rec.seq)) >= 2 * params["window_size"]]

    if not valid_sequences:
        return {"best": best_primers, "other": other_options}

    for name, seq in valid_sequences:
        fw_region = seq[:params["window_size"]]
        rev_region = seq[-params["window_size"]:]

        fw_candidates = find_candidate_primers(fw_region, params["k_range"], params["gc_min"], params["gc_max"])
        # --- CORRECTED REVERSE PRIMER GENERATION ---
        rev_template_candidates = find_candidate_primers(rev_region, params["k_range"], params["gc_min"], params["gc_max"])
        rev_candidates = [reverse_complement(p) for p in rev_template_candidates]

        if not fw_candidates or not rev_candidates:
            best_primers[name] = None
            other_options[name] = []
            continue

        fw_checks = {fw: check_cached_primer_hits(fw, sequences_tuple) for fw in fw_candidates}
        rev_checks = {rv: check_cached_primer_hits(rv, sequences_tuple) for rv in rev_candidates}

        valid_fw = fw_candidates
        valid_rev = rev_candidates
        if mode == "specificity":
            isoform_id = extract_isoform_id(name)
            valid_fw = [
                fw for fw, check in fw_checks.items()
                if isoform_id in check["isoforms"] and check["total"] <= 2
            ]
            valid_rev = [
                rv for rv, check in rev_checks.items()
                if isoform_id in check["isoforms"] and check["total"] <= 2
            ]

        print(f"🔍 [{name}] Found {len(valid_fw)} FW primers, {len(valid_rev)} REV primers after specificity filtering.")

        best_pair = None
        other_pairs = []

        for fw, rv in product(valid_fw, valid_rev):
            primer_pair = {"forward": fw, "reverse": rv, "product_size": None} # Size not applicable for RACE
            if mode == "specificity":
                best_pair = primer_pair
                break
            else:
                if len(other_pairs) < 10:
                    other_pairs.append(primer_pair)
                if not best_pair:
                    best_pair = primer_pair

        # --- REMOVED LOGIC FLAW ---
        # The original code had a 'break' here that prematurely exited the main loop.

        best_primers[name] = best_pair
        other_options[name] = other_pairs

    return {"best": best_primers, "other": other_options}

def find_coverage_primers(sequences, params_qpcr, params_race, primer_type):
    """
    Finds a minimal set of primer pairs to amplify all isoforms (Set Cover Problem).
    This uses a greedy algorithm approach for efficiency.
    """
    all_isoform_names = {rec.id for rec in sequences}
    
    params = params_qpcr if primer_type == "qPCR" else params_race
    
    def get_regions(seq_str):
        if primer_type == "qPCR":
            return (seq_str[params["fw_start"] - 1:params["fw_end"]],
                    seq_str[params["rev_start"] - 1:params["rev_end"]])
        else: # RACE
            return (seq_str[:params["window_size"]], seq_str[-params["window_size"]:])

    # Pre-compute all valid k-mers and the isoforms they hit.
    fw_kmer_hits = {}
    rv_kmer_hits = {}
    for rec in sequences:
        seq_str = str(rec.seq)
        fw_region, rv_region = get_regions(seq_str)
        
        for kmer in extract_kmers(fw_region, params["k_range"]):
            if params["gc_min"] <= calculate_gc_content(kmer) <= params["gc_max"]:
                fw_kmer_hits.setdefault(kmer, set()).add(rec.id)

        for kmer in extract_kmers(rv_region, params["k_range"]):
            if params["gc_min"] <= calculate_gc_content(kmer) <= params["gc_max"]:
                rv_kmer_hits.setdefault(kmer, set()).add(rec.id)

    # Use a greedy algorithm to find the best primer pairs.
    selected_pairs = []
    pending_isoforms = all_isoform_names.copy()
    
    # Create all potential primer pair combinations and their coverage.
    combinations = []
    for fw_kmer, fw_hits in fw_kmer_hits.items():
        for rv_kmer, rv_hits in rv_kmer_hits.items():
            # The coverage of a pair is the intersection of isoforms hit by each primer.
            common_hits = fw_hits.intersection(rv_hits)
            if common_hits:
                combinations.append((fw_kmer, rv_kmer, common_hits))

    while pending_isoforms and combinations:
        # Find the pair that covers the most *remaining* isoforms.
        best_pair = max(combinations, key=lambda item: len(item[2].intersection(pending_isoforms)))
        
        fw, rv, covered_by_best = best_pair
        
        # Check if the best pair actually covers any new isoforms.
        newly_covered = covered_by_best.intersection(pending_isoforms)
        if not newly_covered:
            break # No more progress can be made.

        selected_pairs.append({"forward": fw, "reverse": rv, "isoforms": list(newly_covered)})
        pending_isoforms -= newly_covered # Remove covered isoforms from the pending set.
        
        # Remove the chosen pair to avoid re-selecting it.
        combinations.remove(best_pair)

    return {"selected": selected_pairs, "pending": list(pending_isoforms)}

# ==============================================================================
# --- SCRIPT EXECUTION ---
# ==============================================================================

def create_design(fasta_file_path, primer_type, design_mode):
    """
    Main function to execute the primer design pipeline.
    """
    print(f"📂 Reading FASTA file: {fasta_file_path}")
    try:
        sequence_set = list(SeqIO.parse(fasta_file_path, "fasta"))
        if not sequence_set:
            print("❌ Error: No sequences found in the FASTA file.")
            return
        print(f"✅ Loaded {len(sequence_set)} sequences.")
    except FileNotFoundError:
        print(f"❌ Error: FASTA file not found at '{fasta_file_path}'")
        return

    results = None
    if design_mode == "specificity":
        print(f"🧬 Running in SPECIFICITY mode for {primer_type} primers...")
        if primer_type == "qPCR":
            results = design_qpcr_primers(sequence_set, PARAMS_QPCR, design_mode)
        elif primer_type == "RACE":
            results = design_race_primers(sequence_set, PARAMS_RACE, design_mode)
        else:
            print(f"⚠️ Unsupported primer type: {primer_type}")
            return

    elif design_mode == "coverage":
        print(f"🎯 Running in COVERAGE mode for {primer_type} primers...")
        results = find_coverage_primers(sequence_set, PARAMS_QPCR, PARAMS_RACE, primer_type)
    else:
        print(f"⚠️ Unsupported design mode: {design_mode}")
        return

    # --- Process and Save Results ---
    
    df_rows = []
    if design_mode == "specificity":
        for isoform_name, pair in results["best"].items():
            if not pair: continue
            f, r = pair["forward"], pair["reverse"]
            df_rows.append({
                "Isoform": isoform_name,
                "Primer_Forward": f,
                "GC_Forward": calculate_gc_content(f),
                "Tm_Forward": calculate_melting_temp(f),
                "Primer_Reverse": r,
                "GC_Reverse": calculate_gc_content(r),
                "Tm_Reverse": calculate_melting_temp(r),
                "Hairpin_Forward": find_longest_complementary_run(f, reverse_complement(f)),
                "Hairpin_Reverse": find_longest_complementary_run(r, reverse_complement(r)),
                "SelfDimer_Forward": find_longest_complementary_run(f, f),
                "SelfDimer_Reverse": find_longest_complementary_run(r, r),
                "CrossDimer": find_longest_complementary_run(f, reverse_complement(r)),
            })
    elif design_mode == "coverage":
        for item in results["selected"]:
            df_rows.append({
                "Primer_Forward": item["forward"],
                "Primer_Reverse": item["reverse"],
                "Covered_Isoforms": ";".join(item["isoforms"]),
            })

    if df_rows:
        return pd.DataFrame(df_rows)
    else:
        return None


if __name__ == "__main__":
    FASTA_FILE_PATH = "./data/avengers-3.fa"
    PRIMER_TYPE = "RACE"
    DESIGN_MODE = "coverage" 
    create_design(FASTA_FILE_PATH, PRIMER_TYPE, DESIGN_MODE)
