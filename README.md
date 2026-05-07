# pylcp-nist

`pylcp-nist` is a Python toolkit for modeling laser cooling and trapping physics. It provides building blocks for atomic states, laser beams, magnetic fields, Hamiltonians, rate-equation simulations, optical Bloch equation simulations, and trajectory integration.

This repo is research-oriented and notebook-heavy. The examples are the best place to see how the pieces fit together for molasses, magneto-optical traps, optical pumping, magnetic traps, and related laser-cooling workflows.

## What It Includes

- Atomic state and transition helpers.
- Laser beam, polarization, magnetic-field, and grating utilities.
- Hamiltonian construction tools for internal atomic structure.
- Rate-equation and optical Bloch equation solvers.
- Integration tools for force profiles and particle motion.
- Example notebooks for basic systems, MOTs, molasses, magnetic traps, and laser fields.
- A PDF/manual-style document under `doc/`.

## Project Structure

- `atom.py` - atomic species, states, and transition properties.
- `fields.py` - laser beam and magnetic field definitions.
- `hamiltonian.py` and `hamiltonians/` - Hamiltonian construction helpers.
- `rateeq.py` - rate-equation simulation tools.
- `obe.py` - optical Bloch equation simulation tools.
- `integration_tools.py` - helpers for integrating trajectories and dynamics.
- `examples/` - notebooks and scripts showing common simulation workflows.
- `doc/` - LaTeX/PDF documentation and generated figures.

## Setup

This project currently does not include a `setup.py` or package install configuration. To use it locally, make sure Python can import the parent directory that contains this repo.

One simple local workflow is to add the parent GitHub directory to your Python path, then import the package from notebooks or scripts:

```python
import pylcp
```

Typical dependencies include:

- numpy
- scipy
- matplotlib
- numba

Some notebooks may need additional scientific Python packages depending on the example.

## Where To Start

Start with the notebooks in:

- `examples/basics/`
- `examples/MOTs/`
- `examples/molasses/`

Those examples show the package at a practical level without needing to read through the full solver implementation first.

## Notes

The code appears to be an in-progress research fork, so paths, notebooks, and examples may be more useful than a polished package interface. For deeper background, see `doc/pylcp.pdf`.
