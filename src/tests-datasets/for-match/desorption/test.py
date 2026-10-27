from pathlib import Path

from ase import Atoms
from ase.io import read

p = Path(__file__).parent / "otfkmc.traj"
atoms = read(p)
print(atoms)
assert isinstance(atoms, Atoms)
atoms.write(p.with_name("system.xyz"), format="extxyz")
