# scripts/preprocess.py

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from typing import List
import numpy as np

def smiles_to_morgan_fp(smiles: str, radius: int = 2, n_bits: int = 2048) -> np.ndarray:
    """Convertit une chaîne SMILES en vecteur d'empreinte Morgan (ECFP)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(n_bits)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((n_bits,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def transform_smiles_df(df: pd.DataFrame, smiles_col: str = 'SMILES') -> pd.DataFrame:
    """Applique la transformation SMILES -> empreintes moléculaires pour un DataFrame."""
    fps = df[smiles_col].apply(smiles_to_morgan_fp)
    return pd.DataFrame(fps.tolist(), index=df.index)

def load_and_process_data(path: str) -> pd.DataFrame:
    """Charge et transforme les données SMILES en vecteurs."""
    df = pd.read_csv(path)
    return transform_smiles_df(df)
