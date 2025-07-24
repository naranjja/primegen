import streamlit as st
import tempfile
import os
from primer_designer import create_design


def run_primer_design(fasta_path, primer_type, design_mode):
    return create_design(fasta_path, primer_type, design_mode)


st.title("PRIMEGEN")

uploaded_file = st.file_uploader("Choose FASTA file", type=["fa", "fasta"])
primer_type = st.selectbox("Select Primer Type", ["qPCR", "RACE"])
design_mode = st.selectbox("Select Design Mode", ["specificity", "coverage"])

# Initialize session state
if "output_df" not in st.session_state:
    st.session_state.output_df = None
if "fasta_path" not in st.session_state:
    st.session_state.fasta_path = None

if st.button("Run Primer Design"):
    if uploaded_file is None:
        st.error("Please choose a FASTA file first.")
    else:
        # Save uploaded file to a temporary path
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fa") as tmp_file:
            tmp_file.write(uploaded_file.read())
            st.session_state.fasta_path = tmp_file.name

        # Run the logic
        with st.spinner("Running primer design..."):
            output_df = run_primer_design(
                st.session_state.fasta_path, primer_type, design_mode
            )

        if output_df is None:
            st.error("No suitable primers were found with the given parameters.")
        else:
            st.success("Primer design completed!")
            st.session_state.output_df = output_df

# Show output and download button if results are stored
if st.session_state.output_df is not None:
    st.dataframe(st.session_state.output_df)

    csv = st.session_state.output_df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="Download Results as CSV",
        data=csv,
        file_name="primers_output.csv",
        mime="text/csv",
    )

# Optional cleanup after displaying
if st.session_state.fasta_path and os.path.exists(st.session_state.fasta_path):
    os.remove(st.session_state.fasta_path)
    st.session_state.fasta_path = None