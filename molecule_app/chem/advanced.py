import itertools
import math

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms


def _ideal_angle_deg(atom: Chem.Atom) -> float:
    hybridization = atom.GetHybridization()
    if hybridization == Chem.HybridizationType.SP:
        return 180.0
    if hybridization == Chem.HybridizationType.SP2:
        return 120.0
    return 109.5


def _bond_order_factor(bond: Chem.Bond) -> float:
    bond_type = bond.GetBondTypeAsDouble()
    if bond_type >= 3.0:
        return 0.87
    if bond_type >= 2.0:
        return 0.91
    if bond.GetIsAromatic():
        return 0.93
    return 1.0


def calculate_molecular_strain_energy(
    mol,
    k_bond=300.0,
    k_angle=40.0,
):
    """Estimate strain energy from harmonic bond and angle distortions."""
    conf = mol.GetConformer()
    periodic_table = Chem.GetPeriodicTable()

    bond_terms = []
    total_bond = 0.0
    for bond in mol.GetBonds():
        a_idx = bond.GetBeginAtomIdx()
        b_idx = bond.GetEndAtomIdx()
        atom_a = mol.GetAtomWithIdx(a_idx)
        atom_b = mol.GetAtomWithIdx(b_idx)
        r_measured = rdMolTransforms.GetBondLength(conf, a_idx, b_idx)
        r0 = (
            periodic_table.GetRcovalent(atom_a.GetAtomicNum())
            + periodic_table.GetRcovalent(atom_b.GetAtomicNum())
        ) * _bond_order_factor(bond)
        e_bond = 0.5 * k_bond * (r_measured - r0) ** 2
        total_bond += e_bond
        bond_terms.append(
            {
                "bond": f"{a_idx}-{b_idx}",
                "r_measured_A": round(r_measured, 4),
                "r0_A": round(r0, 4),
                "energy": round(e_bond, 6),
            }
        )

    angle_terms = []
    total_angle = 0.0
    for center_idx in range(mol.GetNumAtoms()):
        center_atom = mol.GetAtomWithIdx(center_idx)
        neighbors = [a.GetIdx() for a in center_atom.GetNeighbors()]
        if len(neighbors) < 2:
            continue
        theta0_deg = _ideal_angle_deg(center_atom)
        theta0 = math.radians(theta0_deg)
        for i_idx, k_idx in itertools.combinations(neighbors, 2):
            theta_deg = rdMolTransforms.GetAngleDeg(conf, i_idx, center_idx, k_idx)
            theta = math.radians(theta_deg)
            e_angle = 0.5 * k_angle * (theta - theta0) ** 2
            total_angle += e_angle
            angle_terms.append(
                {
                    "angle": f"{i_idx}-{center_idx}-{k_idx}",
                    "theta_deg": round(theta_deg, 4),
                    "theta0_deg": theta0_deg,
                    "energy": round(e_angle, 6),
                }
            )

    total_energy = total_bond + total_angle
    return {
        "total_energy": round(total_energy, 6),
        "bond_energy": round(total_bond, 6),
        "angle_energy": round(total_angle, 6),
        "bond_terms": bond_terms,
        "angle_terms": angle_terms,
    }


def predict_cleavage_planes(
    structure,
    pressure_vector,
    max_index=2,
    top_n=6,
):
    """Rank likely cleavage planes from normal stress and interplanar spacing."""
    pressure = np.array(pressure_vector, dtype=float)
    norm = np.linalg.norm(pressure)
    if norm == 0.0:
        raise ValueError("Pressure vector must be non-zero.")
    pressure_dir = pressure / norm

    reciprocal = structure.lattice.reciprocal_lattice_crystallographic.matrix
    number_density = len(structure.sites) / structure.lattice.volume
    results = []

    for h in range(-max_index, max_index + 1):
        for k in range(-max_index, max_index + 1):
            for l in range(-max_index, max_index + 1):
                if h == 0 and k == 0 and l == 0:
                    continue

                g = h * reciprocal[0] + k * reciprocal[1] + l * reciprocal[2]
                g_norm = np.linalg.norm(g)
                if g_norm < 1e-12:
                    continue

                n_hat = g / g_norm
                spacing = 1.0 / g_norm
                resolved_normal = float(np.dot(pressure_dir, n_hat))
                normal_stress = abs(resolved_normal)
                planar_density = number_density * spacing
                cleavage_score = normal_stress * spacing / (planar_density + 1e-12)

                results.append(
                    {
                        "plane": f"({h}{k}{l})",
                        "normal_stress_factor": round(float(normal_stress), 6),
                        "d_spacing_A": round(float(spacing), 6),
                        "estimated_planar_density": round(float(planar_density), 6),
                        "cleavage_score": round(float(cleavage_score), 6),
                    }
                )

    results.sort(key=lambda item: item["cleavage_score"], reverse=True)
    return results[:top_n]
