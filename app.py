import streamlit as st
import tempfile
import os
from primer_designer import create_design


def run_primer_design(fasta_path, primer_type, design_mode):
    return create_design(fasta_path, primer_type, design_mode)


# -------------------------------
# Session State Initialization
# -------------------------------
if "output_df1" not in st.session_state:
    st.session_state.output_df1 = None
if "output_df2" not in st.session_state:
    st.session_state.output_df2 = None
if "fasta_path" not in st.session_state:
    st.session_state.fasta_path = None


# -------------------------------
# Clear output when input changes
# -------------------------------
def clear_output():
    st.session_state.output_df1 = None
    st.session_state.output_df2 = None


# -------------------------------
# UI
# -------------------------------
st.title("🧬 PRIMEGEN")

st.markdown("""
PRIMEGEN is a Python-based primer design tool that identifies optimal primer pairs from multi-isoform FASTA sequences for **qPCR** or **RACE** applications. It operates in two distinct modes:

* **Specificity Mode**: Identifies primer pairs that are unique to a single target isoform, ensuring minimal cross-reactivity. The final output for this mode includes a detailed analysis of primer quality, including GC content, melting temperature (Tm), and predicted secondary structures (hairpins, self-dimers, and cross-dimers) to assess biological viability.

* **Coverage Mode**: Addresses the challenge of amplifying all isoforms with a minimum number of primers. It employs a greedy algorithm to solve this "Set Cover" type problem, prioritizing primer pairs that cover the largest number of previously unamplified transcripts in each step.

The pipeline filters all potential primers based on user-defined parameters for GC content, melting temperature, and k-mer length. For qPCR, it also enforces a specific product size range. Its use of cached functions for primer hit detection ensures efficient and scalable analysis, even with large datasets.
""")

uploaded_file = st.file_uploader(
    "📁 Choose FASTA file", type=["fa", "fasta"], on_change=clear_output
)
primer_type = st.selectbox(
    "🧪 Select Primer Type", ["qPCR", "RACE"], on_change=clear_output
)
design_mode = st.selectbox(
    "🎨 Select Design Mode", ["specificity", "coverage"], on_change=clear_output
)

# -------------------------------
# Run Primer Design
# -------------------------------
if st.button("Run Primer Design ⚡"):
    if uploaded_file is None:
        st.error("Please choose a FASTA file first.")
    else:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fa") as tmp_file:
            tmp_file.write(uploaded_file.read())
            st.session_state.fasta_path = tmp_file.name

        with st.spinner("Running primer design..."):
            output_df1, output_df2 = run_primer_design(
                st.session_state.fasta_path, primer_type, design_mode
            )

        if output_df1 is None:
            st.error("No suitable primers were found with the given parameters.")
        else:
            st.success("Primer design completed!")
            st.session_state.output_df1 = output_df1
            st.session_state.output_df2 = output_df2

# -------------------------------
# Display Results
# -------------------------------
if st.session_state.output_df1 is not None:
    st.subheader("Design")
    st.dataframe(st.session_state.output_df1)

    csv1 = st.session_state.output_df1.to_csv(index=False).encode("utf-8")
    st.download_button(
        label=f"Download Design ({st.session_state.output_df1.shape[0]} pairs)",
        data=csv1,
        file_name=f"primers_design_{primer_type}_{design_mode}.csv",
        mime="text/csv",
    )

if st.session_state.output_df2 is not None:
    st.subheader("Others")
    st.dataframe(st.session_state.output_df2)

    csv2 = st.session_state.output_df2.to_csv(index=False).encode("utf-8")
    st.download_button(
        label=f"Download Others ({st.session_state.output_df2.shape[0]} pairs)",
        data=csv2,
        file_name=f"primers_others_{primer_type}_{design_mode}.csv",
        mime="text/csv",
    )

# -------------------------------
# Clean Up Temp File
# -------------------------------
if st.session_state.fasta_path and os.path.exists(st.session_state.fasta_path):
    os.remove(st.session_state.fasta_path)
    st.session_state.fasta_path = None