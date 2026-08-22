import streamlit as st

from chem.parser import parse_molecule
from chem.geometry import compute_bond_angles, generate_3d_conformer
from chem.crystal import crystal_summary, load_crystal_from_cif
from chem.visualization import molecule_view_html
from chem.advanced import calculate_molecular_strain_energy, predict_cleavage_planes

st.set_page_config(
    page_title="3D Chemistry Math Models",
    page_icon="🧪",
    layout="wide",
)

st.markdown(
    """
<style>
    .stApp {
        background: radial-gradient(circle at top right, #1f2947 0%, #0f172a 45%, #0a0f1f 100%);
        color: #e2e8f0;
    }
    h1, h2, h3, h4, h5, h6, p, label, span {
        color: #e2e8f0 !important;
    }
    .hero-card {
        border: 1px solid rgba(148, 163, 184, 0.25);
        border-radius: 16px;
        padding: 1.25rem 1.5rem;
        background: linear-gradient(135deg, rgba(14, 165, 233, 0.22), rgba(99, 102, 241, 0.22));
        backdrop-filter: blur(3px);
        margin-bottom: 1rem;
    }
    .feature-grid {
        display: grid;
        grid-template-columns: repeat(3, minmax(0, 1fr));
        gap: 0.75rem;
        margin-top: 0.85rem;
    }
    .feature-card {
        border: 1px solid rgba(148, 163, 184, 0.22);
        border-radius: 12px;
        padding: 0.85rem 0.95rem;
        background: rgba(15, 23, 42, 0.65);
    }
    .panel-card {
        border: 1px solid rgba(148, 163, 184, 0.25);
        border-radius: 14px;
        padding: 1rem 1.1rem;
        background: rgba(15, 23, 42, 0.72);
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 0.5rem;
        background: rgba(15, 23, 42, 0.7);
        border-radius: 10px;
        padding: 0.25rem;
    }
    .stTabs [data-baseweb="tab"] {
        border-radius: 8px;
        background: transparent;
    }
    .stTabs [aria-selected="true"] {
        background: rgba(14, 165, 233, 0.2) !important;
        border: 1px solid rgba(125, 211, 252, 0.35) !important;
    }
</style>
""",
    unsafe_allow_html=True,
)

st.markdown(
    """
<div class="hero-card">
  <h1 style="margin:0;">🧬 3D Chemistry Math Models</h1>
  <p style="margin:0.45rem 0 0 0; opacity:0.95;">
    Turn chemistry inputs into interactive 3D molecular geometry, strain-energy estimates, and crystal cleavage predictions.
  </p>
  <div class="feature-grid">
    <div class="feature-card"><b>Molecules</b><br>Build 3D conformers from SMILES/InChI and inspect bond angles.</div>
    <div class="feature-card"><b>Strain energy</b><br>Estimate how much distortion energy is stored in bent or stretched bonds.</div>
    <div class="feature-card"><b>Crystal cleavage</b><br>Rank likely fracture planes under an applied pressure direction.</div>
  </div>
</div>
""",
unsafe_allow_html=True,
)

with st.sidebar:
    st.markdown("### ⚙️ Workspace Guide")
    st.markdown(
        "- Use **SMILES/InChI** for molecular geometry\n"
        "- Use **CIF** for crystal lattice analysis\n"
        "- Advanced tools provide fast numerical estimates for exploration"
    )
    st.info("Tip: start with `CCO` in Molecule and a small CIF sample in Crystal.")

tab1, tab2, tab3 = st.tabs(
    [
        "🧪 Molecule",
        "🧱 Crystal",
        "📐 Advanced Angle",
    ]
)

with tab1:
    c1, c2 = st.columns([2, 1])
    with c1:
        user_input = st.text_input("Enter SMILES or InChI", "CCO")
    with c2:
        st.markdown(
            "<div class='panel-card'><b>Input hints</b><br>Try: CCO, c1ccccc1, InChI strings.</div>",
            unsafe_allow_html=True,
        )

    if st.button("Compute Molecule", type="primary"):
        try:
            mol = parse_molecule(user_input)
            mol = generate_3d_conformer(mol)
            angles = compute_bond_angles(mol)

            m1, m2, m3 = st.columns(3)
            m1.metric("Atoms", mol.GetNumAtoms())
            m2.metric("Bonds", mol.GetNumBonds())
            m3.metric("Angle Samples", len(angles))

            st.subheader("3D Molecule View")
            st.components.v1.html(molecule_view_html(mol), height=520)

            st.subheader("Bond Angles (Sample)")
            st.dataframe(
                [{"i": i, "j": j, "k": k, "angle_deg": a} for i, j, k, a in angles],
                use_container_width=True,
            )
        except Exception as exc:
            st.error(str(exc))

with tab2:
    st.markdown(
        "<div class='panel-card'><b>Crystal mode</b><br>Paste a valid CIF and get lattice constants instantly.</div>",
        unsafe_allow_html=True,
    )
    cif_text = st.text_area("Paste CIF content", height=220)
    if st.button("Analyze Crystal", type="primary"):
        try:
            structure = load_crystal_from_cif(cif_text)
            info = crystal_summary(structure)
            st.subheader("Crystal Lattice Summary")
            s1, s2, s3, s4 = st.columns(4)
            s1.metric("Formula", info["formula"])
            s2.metric("Sites", info["sites"])
            s3.metric("a (A)", info["a"])
            s4.metric("b (A)", info["b"])
            s5, s6, s7 = st.columns(3)
            s5.metric("c (A)", info["c"])
            s6.metric("alpha (deg)", info["alpha"])
            s7.metric("beta / gamma (deg)", f"{info['beta']} / {info['gamma']}")
            with st.expander("Raw summary"):
                st.json(info)
        except Exception as exc:
            st.error(str(exc))

with tab3:
    st.subheader("Custom Strain-Energy + Cleavage-Prediction Algorithms")
    st.caption(
        "These are custom NumPy-based estimators for advanced exploration, "
        "not replacement for full DFT/MD workflows."
    )

    left, right = st.columns(2)
    with left:
        st.markdown("#### Molecular Strain Energy")
        strain_input = st.text_input("SMILES/InChI for strain analysis", "CC1=CC=CC=C1")
    with right:
        st.markdown(
            "<div class='panel-card'><b>What this estimates</b><br>Relative bond/angle distortion energy from current 3D geometry.</div>",
            unsafe_allow_html=True,
        )

    if st.button("Calculate Strain Energy", type="primary"):
        try:
            strain_mol = parse_molecule(strain_input)
            strain_mol = generate_3d_conformer(strain_mol)
            strain_data = calculate_molecular_strain_energy(strain_mol)
            e1, e2, e3 = st.columns(3)
            e1.metric("Total Strain", strain_data["total_energy"])
            e2.metric("Bond Component", strain_data["bond_energy"])
            e3.metric("Angle Component", strain_data["angle_energy"])
            st.markdown("Top bond strain terms")
            st.dataframe(
                sorted(
                    strain_data["bond_terms"],
                    key=lambda term: term["energy"],
                    reverse=True,
                )[:12],
                use_container_width=True,
            )
            st.markdown("Top angle strain terms")
            st.dataframe(
                sorted(
                    strain_data["angle_terms"],
                    key=lambda term: term["energy"],
                    reverse=True,
                )[:12],
                use_container_width=True,
            )
        except ValueError as exc:
            st.error(str(exc))

    st.markdown("#### Crystal Cleavage Under Pressure")
    cleavage_cif = st.text_area(
        "Paste CIF for cleavage prediction",
        key="advanced_cif",
    )
    p_col1, p_col2, p_col3 = st.columns(3)
    with p_col1:
        px = st.number_input("Pressure x", value=1.0)
    with p_col2:
        py = st.number_input("Pressure y", value=0.0)
    with p_col3:
        pz = st.number_input("Pressure z", value=0.0)

    if st.button("Predict Cleavage Planes", type="primary"):
        try:
            cleavage_structure = load_crystal_from_cif(cleavage_cif)
            planes = predict_cleavage_planes(
                cleavage_structure,
                pressure_vector=[px, py, pz],
                max_index=2,
                top_n=8,
            )
            st.success("Predicted likely cleavage planes ranked by custom score.")
            st.dataframe(planes, use_container_width=True)
        except ValueError as exc:
            st.error(str(exc))
