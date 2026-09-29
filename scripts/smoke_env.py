"""Import smoke test for the project environment."""
import botorch
import numpy
import pandas
import sklearn
import torch
from rdkit import Chem, rdBase

mol = Chem.MolFromSmiles("c1ccccc1")
assert mol is not None
assert Chem.MolToSmiles(mol, canonical=True) == "c1ccccc1"

print("RDKit:", rdBase.rdkitVersion)
print("NumPy:", numpy.__version__)
print("Pandas:", pandas.__version__)
print("sklearn:", sklearn.__version__)
print("Torch:", torch.__version__)
print("BoTorch:", botorch.__version__)
print("OK")
