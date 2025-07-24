# PRIMEGEN - Python version of the R script for primer design

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import pandas as pd
from openpyxl import Workbook
import re
import os
from itertools import product

# -------------------------------
# Configuration
# -------------------------------

fasta_path = "./data/avengers-3.fa"
primer_type = "RACE"  # Options: "qPCR" or "RACE"
mode_text = "specificity"  # Options: "specificity" or "coverage"

params_qpcr = {
    "fw_start": 50, "fw_end": 100,
    "rev_start": 150, "rev_end": 200,
    "k_range": range(18, 26),
    "gc_min": 50, "gc_max": 70,
    "prod_min": 80, "prod_max": 150,
    "tm_min": 57, "tm_max": 62
}

params_race = {
    "window_size": 200,
    "k_range": range(23, 29),
    "gc_min": 50, "gc_max": 70,
    "tm_min": 57, "tm_max": 62
}

# -------------------------------
# Utility Functions
# -------------------------------

def calcular_gc(seq):
    gc = sum(1 for b in seq if b in "GC")
    return round(gc / len(seq) * 100, 2) if seq else 0

def extract_kmers(seq, k_range):
    return [seq[i:i+k] for k in k_range for i in range(len(seq) - k + 1)] if seq else []

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def longest_complementary(s, t):
    m, n = len(s), len(t)
    dp = [[0]*(n+1) for _ in range(m+1)]
    max_len = 0
    for i in range(m):
        for j in range(n):
            if is_complement(s[i], t[j]):
                dp[i+1][j+1] = dp[i][j] + 1
                max_len = max(max_len, dp[i+1][j+1])
    return max_len

def is_complement(a, b):
    return (a == "A" and b == "T") or (a == "T" and b == "A") or \
           (a == "C" and b == "G") or (a == "G" and b == "C")

def calcular_tm(seq, Na=50e-3):
    try:
        return round(mt.Tm_NN(seq, Na=Na), 1)
    except:
        return 0.0

def extraer_iso(nombre):
    m = re.search(r"-([0-9]+)", nombre)
    return m.group(1) if m else ""

def verificar_primer(primer, seqs):
    hits = [rec.id for rec in seqs if primer in str(rec.seq)]
    return {"total": len(hits), "isoformas": ",".join(hits)}

def obtener_candidatos(region, k_range, gc_min, gc_max):
    kmers = set(extract_kmers(region, k_range))
    return [k for k in kmers if gc_min <= calcular_gc(k) <= gc_max]

# -------------------------------
# Read input FASTA
# -------------------------------

print("📂 Leyendo archivo FASTA...")
seq_set = list(SeqIO.parse(fasta_path, "fasta"))
print(f"Se han leído {len(seq_set)} secuencias.")

# -------------------------------
# Design logic
# -------------------------------

def generar_qPCR(seqs, params, mode):
    resb = {}
    reso = {}
    for rec in seqs:
        name = rec.id
        seq = str(rec.seq)
        if len(seq) < params["rev_end"]:
            continue

        fwR = seq[params["fw_start"] - 1:params["fw_end"]]
        rvR = seq[params["rev_start"] - 1:params["rev_end"]]

        fw_cands = obtener_candidatos(fwR, params["k_range"], params["gc_min"], params["gc_max"])
        rv_cands = obtener_candidatos(rvR, params["k_range"], params["gc_min"], params["gc_max"])

        best = fallback = None
        others = []

        for fw in fw_cands:
            fchk = verificar_primer(fw, seqs)
            for rv in rv_cands:
                size = params["rev_end"] - params["fw_start"] + len(rv)
                if not (params["prod_min"] <= size <= params["prod_max"]):
                    continue
                rchk = verificar_primer(rv, seqs)
                primer = {"forward": fw, "reverse": rv, "product_size": size}

                valid = (
                    fchk["total"] == 1 and extraer_iso(name) in fchk["isoformas"] and
                    rchk["total"] == 1 and extraer_iso(name) in rchk["isoformas"]
                )
                if mode == "specificity" and valid:
                    best = primer
                    break
                if not fallback:
                    fallback = primer
                if len(others) < 10:
                    others.append(primer)
            if best:
                break

        resb[name] = best or fallback
        reso[name] = others
    return {"best": resb, "other": reso}

def generar_RACE(seqs, params, mode):
    resb = {}
    reso = {}
    for rec in seqs:
        name = rec.id
        seq = str(rec.seq)
        if len(seq) < 2 * params["window_size"]:
            continue

        fwR = seq[:params["window_size"]]
        rvR = seq[-params["window_size"]:]

        fw_cands = obtener_candidatos(fwR, params["k_range"], params["gc_min"], params["gc_max"])
        rv_cands = obtener_candidatos(rvR, params["k_range"], params["gc_min"], params["gc_max"])

        best = fallback = None
        others = []

        for fw in fw_cands:
            fchk = verificar_primer(fw, seqs)
            for rv in rv_cands:
                rchk = verificar_primer(rv, seqs)
                primer = {"forward": fw, "reverse": rv, "product_size": None}

                valid = (
                    fchk["total"] == 1 and extraer_iso(name) in fchk["isoformas"] and
                    rchk["total"] == 1 and extraer_iso(name) in rchk["isoformas"]
                )
                if mode == "specificity" and valid:
                    best = primer
                    break
                if not fallback:
                    fallback = primer
                if len(others) < 10:
                    others.append(primer)
            if best:
                break

        resb[name] = best or fallback
        reso[name] = others
    return {"best": resb, "other": reso}

def generar_coverage(seqs, params_qpcr, params_race, primer_type):
    iso_names = [rec.id for rec in seqs]
    if primer_type == "qPCR":
        get_region = lambda s: (
            extract_kmers(s[params_qpcr["fw_start"] - 1:params_qpcr["fw_end"]], params_qpcr["k_range"]),
            extract_kmers(s[params_qpcr["rev_start"] - 1:params_qpcr["rev_end"]], params_qpcr["k_range"])
        )
        gc_min, gc_max = params_qpcr["gc_min"], params_qpcr["gc_max"]
    else:
        get_region = lambda s: (
            extract_kmers(s[:params_race["window_size"]], params_race["k_range"]),
            extract_kmers(s[-params_race["window_size"]:], params_race["k_range"])
        )
        gc_min, gc_max = params_race["gc_min"], params_race["gc_max"]

    fw_pool, rv_pool = set(), set()
    for rec in seqs:
        s = str(rec.seq)
        fw, rv = get_region(s)
        fw_pool.update(fw)
        rv_pool.update(rv)

    fw_kmers = [k for k in fw_pool if gc_min <= calcular_gc(k) <= gc_max]
    rv_kmers = [k for k in rv_pool if gc_min <= calcular_gc(k) <= gc_max]

    from itertools import product
    candidates = list(product(fw_kmers, rv_kmers))

    def hits(primer):
        return [rec.id for rec in seqs if primer in str(rec.seq)]

    fw_hits = {k: hits(k) for k in fw_kmers}
    rv_hits = {k: hits(k) for k in rv_kmers}

    selected = []
    pending = set(iso_names)
    for fw, rv in candidates:
        common = set(fw_hits[fw]) & set(rv_hits[rv])
        cover = common & pending
        if cover:
            selected.append({"forward": fw, "reverse": rv, "isoforms": list(cover)})
            pending -= cover
            if not pending:
                break

    return {"selected": selected, "pending": list(pending)}

if primer_type == "qPCR":
    if mode_text == "specificity":
        resultados = generar_qPCR(seq_set, params_qpcr, mode_text)
    else:
        resultados = generar_coverage(seq_set, params_qpcr, params_race, primer_type)
else:
    if mode_text == "specificity":
        resultados = generar_RACE(seq_set, params_race, mode_text)
    else:
        resultados = generar_coverage(seq_set, params_qpcr, params_race, primer_type)

# -------------------------------
# Export results to Excel
# -------------------------------

wb = Workbook()
ws_best = wb.active
ws_best.title = "Best_Primers"

headers = ["Isoforma", "Primer_Forward", "Primer_Reverse", "GC_Forward", "GC_Reverse",
           "Tm_Forward", "Tm_Reverse", "Hairpin_Forward", "Hairpin_Reverse",
           "SelfDimer_Forward", "SelfDimer_Reverse", "CrossDimer"]
ws_best.append(headers)

# Handle both specificity and coverage mode outputs
if mode_text == "specificity":
    for nm, p in resultados["best"].items():
        if not p: continue
        f, r = p["forward"], p["reverse"]
        row = [
            nm, f, r,
            calcular_gc(f), calcular_gc(r),
            calcular_tm(f), calcular_tm(r),
            longest_complementary(f, reverse_complement(f)),
            longest_complementary(r, reverse_complement(r)),
            longest_complementary(f, f),
            longest_complementary(r, r),
            longest_complementary(f, reverse_complement(r))
        ]
        ws_best.append(row)
else:
    ws_cov = wb.create_sheet("Coverage_Primers")
    ws_cov.append(["Primer_Forward", "Primer_Reverse", "Isoformas"])
    for item in resultados["selected"]:
        ws_cov.append([item["forward"], item["reverse"], ";".join(item["isoforms"])])

# Create output directory if it doesn't exist
os.makedirs("./output", exist_ok=True)
output_path = "./output/results.xlsx"
wb.save(output_path)
print(f"✅ Resultados guardados en: {output_path}")