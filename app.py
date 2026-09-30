import streamlit as st
import streamlit.components.v1 as components
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, MeltingTemp as mt
import pandas as pd

# Page Configuration
st.set_page_config(page_title="BioMedical Sequence & Diagnostic Analyzer", page_icon="🧬", layout="wide")

st.title("🧬 BioMedical Sequence & Diagnostic Analyzer (v3.1)")
st.write("An engineering-focused bioinformatics platform featuring direct NCBI genomic retrieval, 3D macromolecular visualization, and clinical diagnostic workflows.")

# Initialize session state for sequence sharing across tabs
if 'shared_sequence' not in st.session_state:
    st.session_state['shared_sequence'] = "ATGCGATCGATCGATCGATCGATCGATCGATCTAG"

# Application Navigation Tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🌐 NCBI Direct Fetcher",
    "📊 Basic Analysis & Translation", 
    "🔬 Diagnostic PCR Primer Design", 
    "🩺 Clinical Mutation Detector",
    "🧊 3D Protein Structure Viewer"
])

# ---------------------------------------------------------
# TAB 1: NCBI Database Integration (Entrez API)
# ---------------------------------------------------------
with tab1:
    st.header("🌐 NCBI Entrez Database Direct Fetcher")
    st.write("Fetch nucleotide sequences directly from the National Center for Biotechnology Information (NCBI) database using Accession IDs.")
    
    col_entrez_1, col_entrez_2 = st.columns([2, 1])
    with col_entrez_1:
        accession_id = st.text_input("Enter NCBI Accession Number / ID:", "NM_000207", help="Examples: NM_000207 (Human Insulin), NC_045512 (SARS-CoV-2 genome)").strip()
    with col_entrez_2:
        user_email = st.text_input("Enter Email Address (NCBI Requirement):", "researcher@example.com").strip()
        
    fetch_button = st.button("🔍 Fetch Sequence from NCBI", use_container_width=True)
    
    if fetch_button:
        if not accession_id:
            st.error("⚠️ Please enter a valid NCBI Accession Number.")
        elif not user_email or "@" not in user_email:
            st.error("⚠️ Please enter a valid email address required by NCBI Entrez policy.")
        else:
            try:
                with st.spinner(f"Fetching record '{accession_id}' from NCBI Nucleotide database..."):
                    Entrez.email = user_email
                    handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
                    record = SeqIO.read(handle, "fasta")
                    handle.close()
                    
                    st.success(f"✅ Successfully fetched: {record.id}")
                    st.subheader("📋 Sequence Metadata")
                    st.write(f"**Description:** {record.description}")
                    
                    m_col1, m_col2 = st.columns(2)
                    m_col1.metric("Sequence Length (bp)", len(record.seq))
                    m_col2.metric("GC-Content", f"{gc_fraction(record.seq)*100:.2f}%")
                    
                    fetched_str = str(record.seq).upper()
                    st.subheader("🧬 FASTA Nucleotide Sequence")
                    st.text_area("Fetched DNA Sequence:", fetched_str, height=180, key="fetched_seq_box")
                    
                    st.session_state['shared_sequence'] = fetched_str
                    st.info("💡 This sequence has been automatically transferred to the 'Basic Analysis & Translation' tab!")
                    
            except Exception as e:
                st.error(f"❌ Failed to retrieve record from NCBI. Error: {str(e)}")

# ---------------------------------------------------------
# TAB 2: Basic Analysis & Translation
# ---------------------------------------------------------
with tab2:
    st.header("📊 Sequence Analysis & Central Dogma")
    sequence_input = st.text_area(
        "Enter DNA sequence (or edit fetched NCBI sequence):", 
        value=st.session_state['shared_sequence'], 
        height=150, 
        key="t2_seq"
    )
    cleaned_seq = "".join(sequence_input.split()).upper()

    if cleaned_seq:
        if set(cleaned_seq).issubset(set("ATCGN")):
            dna_seq = Seq(cleaned_seq)
            col1, col2 = st.columns(2)
            col1.metric("Sequence Length (bp)", len(dna_seq))
            col2.metric("GC-Content", f"{gc_fraction(dna_seq)*100:.2f}%")
            
            st.subheader("Transcription & Translation")
            st.write("**RNA Sequence:**")
            st.code(str(dna_seq.transcribe()), language="text")
            
            protein = dna_seq.translate()
            st.write("**Protein Sequence (Amino Acids):**")
            st.code(str(protein), language="text")
            
            st.subheader("📈 Amino Acid Frequency Distribution")
            aa_counts = {aa: protein.count(aa) for aa in set(protein) if aa != "*"}
            if aa_counts:
                df = pd.DataFrame(list(aa_counts.items()), columns=['Amino Acid', 'Frequency']).set_index('Amino Acid')
                st.bar_chart(df)
        else:
            st.error("⚠️ Invalid DNA sequence! Please use A, T, C, G bases only.")

# ---------------------------------------------------------
# TAB 3: Diagnostic PCR Primer Design
# ---------------------------------------------------------
with tab3:
    st.header("🔬 Diagnostic PCR Primer Thermodynamics")
    st.write("Analyze primer suitability for thermal cyclers and diagnostic qPCR devices.")
    
    primer_input = st.text_input("Enter Diagnostic Primer Sequence (18-30 bp):", "ATGCGATCGATCGATCGATC", key="t3_primer")
    cleaned_primer = "".join(primer_input.split()).upper()
    
    if cleaned_primer:
        if set(cleaned_primer).issubset(set("ATCG")):
            primer_seq = Seq(cleaned_primer)
            
            tm_val = mt.Tm_NN(primer_seq)
            gc_val = gc_fraction(primer_seq) * 100
            
            c1, c2, c3 = st.columns(3)
            c1.metric("Primer Length", f"{len(primer_seq)} bp")
            c2.metric("GC Ratio", f"{gc_val:.1f}%")
            c3.metric("Melting Temp ($T_m$)", f"{tm_val:.2f} °C")
            
            st.subheader("Diagnostic Suitability Verdict:")
            if 18 <= len(primer_seq) <= 30 and 55 <= tm_val <= 65 and 40 <= gc_val <= 60:
                st.success("✅ **Optimal Primer:** Excellent parameters for standard diagnostic PCR assays.")
            else:
                st.warning("⚠️ **Sub-optimal Primer:** Parameters fall outside ideal qPCR device limits ($T_m$: 55-65°C, GC: 40-60%, Length: 18-30 bp).")

# ---------------------------------------------------------
# TAB 4: Clinical Mutation Detector
# ---------------------------------------------------------
with tab4:
    st.header("🩺 Clinical Mutation & Variant Detector")
    st.write("Compare a patient sequence against a reference sequence to detect single nucleotide mutations (SNPs).")
    
    ref_seq = st.text_input("Reference DNA Sequence (Normal):", "ATGCGATCGATCGATC", key="ref_seq").strip().upper()
    pat_seq = st.text_input("Patient DNA Sequence (Sample):", "ATGCGACCGATCGATC", key="pat_seq").strip().upper()
    
    if ref_seq and pat_seq:
        if len(ref_seq) != len(pat_seq):
            st.warning("⚠️ Sequences must be of equal length for direct alignment comparison.")
        else:
            mutations = []
            for i in range(len(ref_seq)):
                if ref_seq[i] != pat_seq[i]:
                    mutations.append((i + 1, ref_seq[i], pat_seq[i]))
            
            if mutations:
                st.error(f"🚨 Detected {len(mutations)} mutation(s) in the patient sample:")
                mut_df = pd.DataFrame(mutations, columns=["Position (bp)", "Reference Base", "Patient Base"])
                st.table(mut_df)
            else:
                st.success("✅ No point mutations detected! Patient sequence matches reference 100%.")

# ---------------------------------------------------------
# TAB 5: 3D Protein Structure Viewer (Mobile Optimized)
# ---------------------------------------------------------
with tab5:
    st.header("🧊 3D Macromolecular & Protein Structure Viewer")
    st.write("Render 3D protein structures from PDB. *Optimized for mobile GPUs.*")
    
    p_col1, p_col2 = st.columns([2, 1])
    with p_col1:
        pdb_id = st.text_input("Enter 4-character PDB ID:", "1CRN", help="Recommended light structures: 1CRN (Ultra-light), 1TCO, 4HHB").strip().upper()
    with p_col2:
        style_choice = st.selectbox("Style:", ["cartoon", "stick", "sphere", "line"])
        
    spin_toggle = st.checkbox("Enable 3D Auto-Rotation (Uses extra GPU)", value=False)
    render_btn = st.button("🚀 Render 3D Model", use_container_width=True)
    
    if render_btn and pdb_id:
        if len(pdb_id) == 4:
            spin_code = "viewer.spin(true);" if spin_toggle else ""
            html_3d_code = f"""
            <div id="container" style="width: 100%; height: 380px; position: relative;"></div>
            <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
            <script>
                let element = document.getElementById('container');
                let config = {{ backgroundColor: 'white' }};
                let viewer = $3Dmol.createViewer(element, config);
                $3Dmol.download('pdb:{pdb_id}', viewer, {{}}, function() {{
                    viewer.setStyle({{}}, {{{style_choice}: {{color: 'spectrum'}}}});
                    {spin_code}
                    viewer.render();
                }});
            </script>
            """
            components.html(html_3d_code, height=400)
            st.caption(f"📍 Rendered PDB Entry: **{pdb_id}** from RCSB Protein Data Bank.")
        else:
            st.warning("⚠️ PDB IDs must be exactly 4 characters.")
