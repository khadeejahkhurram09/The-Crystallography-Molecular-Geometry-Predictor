import re

import streamlit as st

from chem.advanced import calculate_molecular_strain_energy, predict_cleavage_planes
from chem.crystal import crystal_summary, load_crystal_from_cif
from chem.geometry import compute_bond_angles, generate_3d_conformer
from chem.parser import parse_molecule
from chem.visualization import molecule_view_html


MOLECULE_CHARS = re.compile(r"^[A-Za-z0-9@+\-\[\]\(\)=#\\/%.,:]+$")
FORMULA_ONLY = re.compile(r"^[A-Za-z0-9]+$")


def _normalize_molecule_input(raw_input: str) -> str:
    text = raw_input.strip()
    if not text:
        raise ValueError("Please enter a SMILES or InChI string.")
    if text.lower().startswith("inchi="):
        text = "InChI=" + text.split("=", 1)[1]
    return text


def _molecule_input_message(raw_input: str) -> str:
    text = raw_input.strip()
    if not text:
        return "Please enter a SMILES or InChI string."
    if FORMULA_ONLY.fullmatch(text):
        if text.lower() == "h2o":
            return (
                "That looks like the formula H2O. This app needs a structural string, "
                "so try `O` or `InChI=1S/H2O/h1H2` instead."
            )
        return (
            "That looks like a chemical formula, not a structural string. "
            "Use SMILES or InChI so the app can build a 3D model."
        )
    if not MOLECULE_CHARS.fullmatch(text):
        return "Your molecule input contains unsupported characters. Use only SMILES or InChI syntax."
    return "The molecule string could not be parsed. Check for typos and try a valid SMILES or InChI string."


def _run_molecule_analysis(raw_input: str):
    normalized = _normalize_molecule_input(raw_input)
    if FORMULA_ONLY.fullmatch(normalized):
        raise ValueError(_molecule_input_message(normalized))
    if not MOLECULE_CHARS.fullmatch(normalized):
        raise ValueError(_molecule_input_message(normalized))

    try:
        mol = parse_molecule(normalized)
        mol = generate_3d_conformer(mol)
        angles = compute_bond_angles(mol)
    except ValueError as exc:
        message = str(exc)
        if "Invalid molecule input" in message:
            raise ValueError(_molecule_input_message(normalized)) from exc
        if "3D embedding failed" in message:
            raise ValueError(
                "The molecule was recognized, but 3D embedding failed. Try a simpler or different valid SMILES/InChI string."
            ) from exc
        raise

    return mol, angles


def _run_crystal_analysis(cif_text: str):
    text = cif_text.strip()
    if not text:
        raise ValueError("Please paste CIF content before running crystal analysis.")
    try:
        structure = load_crystal_from_cif(text)
    except Exception as exc:
        raise ValueError("The CIF text could not be parsed. Please paste valid CIF content.") from exc
    return crystal_summary(structure)


def _run_strain_analysis(raw_input: str):
    mol, _ = _run_molecule_analysis(raw_input)
    return calculate_molecular_strain_energy(mol)


def _run_cleavage_analysis(cif_text: str, pressure_vector):
    text = cif_text.strip()
    if not text:
        raise ValueError("Please paste CIF content before running cleavage prediction.")
    try:
        structure = load_crystal_from_cif(text)
    except Exception as exc:
        raise ValueError("The CIF text could not be parsed. Please paste valid CIF content.") from exc
    return predict_cleavage_planes(structure, pressure_vector=pressure_vector, max_index=2, top_n=8)


def _store_state(prefix: str, *, input_value=None, result=None, error=None, extra=None):
    if input_value is not None:
        st.session_state[f"{prefix}_input"] = input_value
    st.session_state[f"{prefix}_result"] = result
    st.session_state[f"{prefix}_error"] = error
    if extra:
        for key, value in extra.items():
            st.session_state[f"{prefix}_{key}"] = value


def _render_molecule_result(result):
    mol = result["mol"]
    angles = result["angles"]

    col1, col2, col3 = st.columns(3)
    col1.metric("Atoms", mol.GetNumAtoms())
    col2.metric("Bonds", mol.GetNumBonds())
    col3.metric("Angle Samples", len(angles))

    st.subheader("3D Molecule View")
    st.components.v1.html(molecule_view_html(mol), height=520)

    st.subheader("Bond Angles")
    if angles:
        st.metric("First Sample Angle", f"{angles[0][3]} deg")
        st.dataframe(
            [{"i": i, "j": j, "k": k, "angle_deg": a} for i, j, k, a in angles],
            use_container_width=True,
        )
    else:
        st.info("No bond angles were found for this structure.")


def _render_crystal_result(info):
    row1 = st.columns(4)
    row1[0].metric("Formula", info["formula"])
    row1[1].metric("Sites", info["sites"])
    row1[2].metric("a (A)", info["a"])
    row1[3].metric("b (A)", info["b"])

    row2 = st.columns(3)
    row2[0].metric("c (A)", info["c"])
    row2[1].metric("alpha (deg)", info["alpha"])
    row2[2].metric("beta / gamma (deg)", f"{info['beta']} / {info['gamma']}")

    with st.expander("Raw summary"):
        st.json(info)


def _render_strain_result(strain_data):
    row1 = st.columns(3)
    row1[0].metric("Total Strain", strain_data["total_energy"])
    row1[1].metric("Bond Component", strain_data["bond_energy"])
    row1[2].metric("Angle Component", strain_data["angle_energy"])

    st.markdown("Top bond strain terms")
    st.dataframe(
        sorted(strain_data["bond_terms"], key=lambda term: term["energy"], reverse=True)[:12],
        use_container_width=True,
    )

    st.markdown("Top angle strain terms")
    st.dataframe(
        sorted(strain_data["angle_terms"], key=lambda term: term["energy"], reverse=True)[:12],
        use_container_width=True,
    )


def _render_cleavage_result(planes):
    if not planes:
        st.info("No candidate cleavage planes were generated.")
        return

    top = planes[0]
    row1 = st.columns(4)
    row1[0].metric("Top Plane", top["plane"])
    row1[1].metric("Top Score", top["cleavage_score"])
    row1[2].metric("d-spacing (A)", top["d_spacing_A"])
    row1[3].metric("Normal Stress", top["normal_stress_factor"])

    st.dataframe(planes, use_container_width=True)


def main():
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
        st.markdown("### Controls")

        with st.form("molecule_form"):
            st.markdown("#### Molecule")
            molecule_input = st.text_input(
                "SMILES or InChI",
                value=st.session_state.get("molecule_input", "CCO"),
                help="Use SMILES or InChI. Example: CCO or InChI=1S/H2O/h1H2",
            )
            molecule_submit = st.form_submit_button("Compute Molecule")

        with st.form("crystal_form"):
            st.markdown("#### Crystal")
            crystal_input = st.text_area(
                "Paste CIF content",
                value=st.session_state.get("crystal_input", ""),
                height=160,
            )
            crystal_submit = st.form_submit_button("Analyze Crystal")

        with st.form("strain_form"):
            st.markdown("#### Strain Energy")
            strain_input = st.text_input(
                "SMILES or InChI for strain analysis",
                value=st.session_state.get("strain_input", "CC1=CC=CC=C1"),
            )
            strain_submit = st.form_submit_button("Calculate Strain Energy")

        with st.form("cleavage_form"):
            st.markdown("#### Cleavage Prediction")
            cleavage_input = st.text_area(
                "Paste CIF for cleavage prediction",
                value=st.session_state.get("cleavage_input", ""),
                height=160,
            )
            px = st.number_input("Pressure x", value=float(st.session_state.get("pressure_x", 1.0)))
            py = st.number_input("Pressure y", value=float(st.session_state.get("pressure_y", 0.0)))
            pz = st.number_input("Pressure z", value=float(st.session_state.get("pressure_z", 0.0)))
            cleavage_submit = st.form_submit_button("Predict Cleavage Planes")

        st.markdown(
            "- Use **SMILES/InChI** for molecular geometry\n"
            "- Use **CIF** for crystal lattice analysis\n"
            "- Advanced tools provide fast numerical estimates for exploration"
        )

    if molecule_submit:
        try:
            mol, angles = _run_molecule_analysis(molecule_input)
            _store_state(
                "molecule",
                input_value=molecule_input,
                result={"mol": mol, "angles": angles},
                error=None,
            )
        except ValueError as exc:
            _store_state("molecule", input_value=molecule_input, result=None, error=str(exc))

    if crystal_submit:
        try:
            info = _run_crystal_analysis(crystal_input)
            _store_state("crystal", input_value=crystal_input, result=info, error=None)
        except ValueError as exc:
            _store_state("crystal", input_value=crystal_input, result=None, error=str(exc))

    if strain_submit:
        try:
            strain_data = _run_strain_analysis(strain_input)
            _store_state("strain", input_value=strain_input, result=strain_data, error=None)
        except ValueError as exc:
            _store_state("strain", input_value=strain_input, result=None, error=str(exc))

    if cleavage_submit:
        try:
            planes = _run_cleavage_analysis(cleavage_input, [px, py, pz])
            _store_state(
                "cleavage",
                input_value=cleavage_input,
                result=planes,
                error=None,
                extra={"pressure_x": px, "pressure_y": py, "pressure_z": pz},
            )
        except ValueError as exc:
            _store_state(
                "cleavage",
                input_value=cleavage_input,
                result=None,
                error=str(exc),
                extra={"pressure_x": px, "pressure_y": py, "pressure_z": pz},
            )

    tab1, tab2, tab3 = st.tabs(["🧪 Molecule", "🧱 Crystal", "📐 Advanced Angle"])

    with tab1:
        st.subheader("Molecule Builder")
        st.caption("Enter a structural string in the sidebar and render the 3D model here.")
        if st.session_state.get("molecule_error"):
            st.error(st.session_state["molecule_error"])
        elif st.session_state.get("molecule_result"):
            _render_molecule_result(st.session_state["molecule_result"])
        else:
            st.info("Use the sidebar to compute a molecule. The 3D model will appear here.")

    with tab2:
        st.subheader("Crystal Lattice Analyzer")
        st.caption("Paste CIF content in the sidebar to inspect lattice geometry.")
        if st.session_state.get("crystal_error"):
            st.error(st.session_state["crystal_error"])
        elif st.session_state.get("crystal_result"):
            _render_crystal_result(st.session_state["crystal_result"])
        else:
            st.info("Use the sidebar to analyze a crystal structure.")

    with tab3:
        st.subheader("Custom Strain-Energy + Cleavage Prediction")
        st.caption("These are fast geometry-based estimators for exploration, not full simulation workflows.")

        if st.session_state.get("strain_error"):
            st.error(st.session_state["strain_error"])
        elif st.session_state.get("strain_result"):
            _render_strain_result(st.session_state["strain_result"])
        else:
            st.info("Use the sidebar to calculate molecular strain energy.")

        st.markdown("#### Crystal Cleavage Under Pressure")
        if st.session_state.get("cleavage_error"):
            st.error(st.session_state["cleavage_error"])
        elif st.session_state.get("cleavage_result"):
            _render_cleavage_result(st.session_state["cleavage_result"])
        else:
            st.info("Use the sidebar to predict likely cleavage planes.")


if __name__ == "__main__":
    main()
