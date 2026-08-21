from pymatgen.core import Structure


def load_crystal_from_cif(cif_text: str):
    """Parse crystal structure from CIF text."""
    return Structure.from_str(cif_text, fmt="cif")


def crystal_summary(structure):
    lat = structure.lattice
    return {
        "a": round(lat.a, 4),
        "b": round(lat.b, 4),
        "c": round(lat.c, 4),
        "alpha": round(lat.alpha, 2),
        "beta": round(lat.beta, 2),
        "gamma": round(lat.gamma, 2),
        "formula": structure.composition.reduced_formula,
        "sites": len(structure.sites),
    }
