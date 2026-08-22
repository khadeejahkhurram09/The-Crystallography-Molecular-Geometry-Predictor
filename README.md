# The Crystallography Molecular Geometry Predictor

An interactive Streamlit application for visualizing molecular geometry and crystal structure in 3D, with custom estimators for strain energy and cleavage behavior.

It converts chemical inputs into interactive models and numerical summaries so users can explore how atoms are arranged, how bonds distort, and how crystals may fracture under pressure.

## Project description

Traditional chemistry diagrams are often static and two-dimensional. This app brings those ideas into a more realistic 3D environment by combining molecular modeling, crystallographic analysis, and lightweight geometry-based heuristics.

It is designed for:

- students learning molecular shape and crystal structure
- researchers and hobbyists exploring structures interactively
- anyone who wants a clearer geometric view of chemistry data

The app is not a replacement for full quantum chemistry or production-grade materials simulation. Instead, it provides a fast, visual, and intuitive way to understand structure-property relationships.

## What this app does

This project helps users explore molecules and crystals in 3D instead of only seeing them as static textbook diagrams.

It supports three main workflows:

- **Molecule builder**: enter a SMILES or InChI string and generate a 3D molecular structure
- **Crystal analysis**: paste CIF data and inspect lattice constants and crystal information
- **Advanced estimation tools**:
  - estimate **molecular strain energy**
  - predict likely **crystal cleavage planes** under pressure

## Why it exists

Chemistry is often taught using simplified 2D drawings. This app shows how the same structure behaves in real 3D space:

- atoms occupy specific positions
- bonds form angles and distances
- molecules can bend, twist, and store strain energy
- crystals can fracture along geometric planes

The goal is to make these ideas easier to visualize and explore.

## Features

### 1) 3D Molecule Viewer

- Accepts **SMILES** or **InChI**
- Adds hydrogens and generates a 3D conformer
- Displays the molecule interactively in 3D
- Samples bond angles from the geometry

### 2) Crystal Lattice Analyzer

- Accepts **CIF** crystal structure text
- Parses the lattice with Pymatgen
- Displays:
  - formula
  - number of sites
  - lattice constants `a`, `b`, `c`
  - lattice angles `alpha`, `beta`, `gamma`

### 3) Molecular Strain Energy

This is a custom numerical estimator that approximates how much distortion energy is stored in a molecule.

It uses:

- bond stretching terms
- bond angle bending terms

The result is a fast heuristic, not a replacement for full quantum chemistry.

### 4) Crystal Cleavage Prediction

This feature estimates which crystal planes are more likely to cleave under an applied pressure direction.

It uses:

- reciprocal lattice geometry
- plane spacing
- pressure direction alignment
- a simple planar-density heuristic

Again, this is an exploratory model rather than a full materials simulation pipeline.

## How it works

### Molecular geometry

1. Parse the input molecule with RDKit
2. Add hydrogens
3. Generate a 3D conformer
4. Optimize the geometry
5. Measure bond angles and display the structure

### Strain energy

The app compares:

- measured bond lengths vs. estimated equilibrium bond lengths
- measured bond angles vs. idealized hybridization angles

It then uses harmonic energy terms to estimate total strain.

### Crystal cleavage

The app loops over Miller indices `(hkl)` and ranks possible planes using a score based on:

- how strongly pressure is aligned with the plane normal
- how widely spaced the planes are
- how dense the plane is

## Algorithmic decisions

The app uses a few intentional mathematical choices rather than arbitrary transforms:

### 1) Reciprocal lattice matrix for crystal planes

For cleavage prediction, the code builds plane normals from the **reciprocal lattice** instead of the direct lattice. That is the correct space for Miller indices `(hkl)` because a plane index represents orientation and spacing, not atomic coordinates directly.

- direct lattice vectors describe the real-space unit cell
- reciprocal lattice vectors describe plane normals and interplanar spacing
- using the reciprocal matrix makes the `(hkl)` search mathematically consistent and easy to rank

### 2) Harmonic energy for strain

Bond stretching and bond bending are modeled with quadratic terms:

- small distortions give small energy changes
- larger distortions rise faster
- the result behaves like a spring approximation, which is a common first-order physical model

This keeps the estimator fast and interpretable.

### 3) Normalized pressure direction

Pressure is converted into a unit vector before scoring cleavage planes. This separates:

- **direction**, which affects which plane is favored
- **magnitude**, which is not needed for the ranking heuristic

That makes the score easier to compare across inputs.

### 4) Heuristic ranking instead of full simulation

The cleavage model uses geometric proxies such as spacing and planar density rather than a full fracture mechanics solver. That choice is intentional because the app is meant to be lightweight, fast, and educational.

## Tech stack

- **Streamlit** for the UI
- **RDKit** for molecular parsing, conformers, and bond geometry
- **Pymatgen** for crystal parsing and lattice data
- **NumPy** for vector and matrix math
- **py3Dmol** for interactive 3D visualization

## Installation

Install dependencies:

```bash
pip install -r requirements.txt
```

## Run the app

```bash
streamlit run app.py
```

## Example inputs

### Molecule

- `CCO`
- `c1ccccc1`
- `InChI=1S/H2O/h1H2`

### Crystal

Paste valid CIF content into the crystal text box.

## Project structure

- `app.py` - Streamlit user interface
- `chem/parser.py` - molecule parsing
- `chem/geometry.py` - 3D conformer generation and bond-angle sampling
- `chem/crystal.py` - CIF parsing and crystal summaries
- `chem/advanced.py` - strain-energy and cleavage heuristics
- `chem/visualization.py` - 3D molecule rendering

## Notes

- The advanced calculations are designed for exploration and education.
- They are not a substitute for DFT, molecular dynamics, or full crystallographic simulation workflows.
- Crystal cleavage prediction is a heuristic ranking model, not a laboratory-grade fracture solver.
