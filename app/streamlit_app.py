# app/streamlit_app.py

import streamlit as st
import pandas as pd
import numpy as np
import joblib
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from PIL import Image

# ----------- Paramètres de base ----------- #
st.set_page_config(page_title="ADMET Predictor", page_icon="💊")

st.title("💊 Prédiction des Propriétés ADMET")
st.markdown("Ce prototype utilise un modèle de machine learning pour prédire 3 propriétés ADMET à partir d'une molécule (SMILES).")

# ----------- Fonctions ----------- #

def smiles_to_morgan_fp(smiles: str, radius: int = 2, n_bits: int = 2048) -> np.ndarray:
    """Convertit une molécule SMILES en empreinte moléculaire binaire (vecteur numpy)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(n_bits)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def draw_molecule(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Draw.MolToImage(mol, size=(300, 300))
    return None

@st.cache_resource
def load_model(path: str):
    return joblib.load(path)

# ----------- Interface utilisateur ----------- #

smiles_input = st.text_input("👉 Entrez une chaîne SMILES :", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")

if smiles_input:
    mol_img = draw_molecule(smiles_input)
    if mol_img:
        st.image(mol_img, caption="Structure moléculaire", use_column_width=False)
    else:
        st.error("⚠️ SMILES invalide. Veuillez réessayer.")

    # ----------- Prédiction ----------- #
    model = load_model("models/admet_model.pkl")
    features = smiles_to_morgan_fp(smiles_input).reshape(1, -1)

    if model:
        prediction = model.predict(features)[0]
        labels = ["Absorption optimale (Y1)", "Distribution optimale (Y2)", "Métabolisme/Excrétion non toxique (Y3)"]
        st.subheader("🧠 Prédictions ADMET :")
        for label, value in zip(labels, prediction):
            st.success(f"{label} : {'✅ Oui' if value == 1 else '❌ Non'}")
    else:
        st.error("Modèle introuvable.")
