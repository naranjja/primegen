import streamlit as st
import pandas as pd
import tempfile
import os

# Placeholder for your primer design function
def run_primer_design(fasta_path, primer_type, design_mode):
    # Replace this with your real logic
    # Simulated output:
    data = {
        "Primer ID": ["primer1", "primer2"],
        "Sequence": ["ATCGTAGCTAGCT", "CGATCGATGCTA"],
        "Tm": [60.1, 59.5],
        "GC%": [52, 48]
    }
    df = pd.DataFrame(data)
    return df

# Streamlit UI
st.title("Primer Design App")

uploaded_file = st.file_uploader("Upload FASTA file", type=["fa", "fasta"])
primer_type = st.selectbox("Select Primer Type", ["qPCR", "RACE"])
design_mode = st.selectbox("Select Design Mode", ["specificity", "coverage"])

if st.button("Run Primer Design"):
    if uploaded_file is None:
        st.error("Please upload a FASTA file first.")
    else:
        # Save uploaded file to a temporary path
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fa") as tmp_file:
            tmp_file.write(uploaded_file.read())
            fasta_path = tmp_file.name

        # Run the logic
        with st.spinner("Running primer design..."):
            output_df = run_primer_design(fasta_path, primer_type, design_mode)

        st.success("Primer design completed!")
        st.dataframe(output_df)

        # Create CSV for download
        csv = output_df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="Download Results as CSV",
            data=csv,
            file_name="primers_output.csv",
            mime="text/csv"
        )

        # Clean up temp file
        os.remove(fasta_path)
