# 🧬 PRIMEGEN
PRIMEGEN is a Python-based primer design tool that identifies optimal primer pairs from multi-isoform FASTA sequences for **qPCR** or **RACE** applications. It operates in two distinct modes:

* **Specificity Mode**: Identifies primer pairs that are unique to a single target isoform, ensuring minimal cross-reactivity. The final output for this mode includes a detailed analysis of primer quality, including GC content, melting temperature (Tm), and predicted secondary structures (hairpins, self-dimers, and cross-dimers) to assess biological viability.

* **Coverage Mode**: Addresses the challenge of amplifying all isoforms with a minimum number of primers. It employs a greedy algorithm to solve this "Set Cover" type problem, prioritizing primer pairs that cover the largest number of previously unamplified transcripts in each step.

The pipeline filters all potential primers based on user-defined parameters for GC content, melting temperature, and k-mer length. For qPCR, it also enforces a specific product size range. Its use of cached functions for primer hit detection ensures efficient and scalable analysis, even with large datasets.

### Run it locally:
- `virtualenv venv`
- `source venv/bin/activate` (UNIX) or `.\venv\Scripts\Activate` (Windows)
- `pip install -r requirements.txt`
- `streamlit run app.py`