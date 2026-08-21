import py3Dmol
from rdkit.Chem import MolToMolBlock


def molecule_view_html(mol):
    """Return an HTML widget for rendering the molecule in 3D."""
    mb = MolToMolBlock(mol)
    view = py3Dmol.view(width=700, height=500)
    view.addModel(mb, "mol")
    view.setStyle({"stick": {}})
    view.addSurface(py3Dmol.VDW, {"opacity": 0.2})
    view.zoomTo()
    return view._make_html()
