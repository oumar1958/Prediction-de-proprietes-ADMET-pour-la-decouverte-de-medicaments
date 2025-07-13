# scripts/train.py

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.multioutput import MultiOutputClassifier
from sklearn.metrics import f1_score
import joblib

from scripts.preprocess import load_and_process_data

def train_model(x_path: str, y_path: str, model_output_path: str):
    X_train_raw = pd.read_csv(x_path)
    y_train = pd.read_csv(y_path).drop(columns=["ID"])
    
    print("Transformation des SMILES en empreintes...")
    X_train = load_and_process_data(x_path)

    print("Entraînement du modèle...")
    base_model = RandomForestClassifier(n_estimators=100, random_state=42)
    model = MultiOutputClassifier(base_model)
    model.fit(X_train, y_train)

    print("Modèle entraîné. Sauvegarde...")
    joblib.dump(model, model_output_path)
    print(f"Modèle sauvegardé sous {model_output_path}")

if __name__ == '__main__':
    train_model('data/raw/X_train.csv', 'data/raw/y_train.csv', 'models/admet_model.pkl')
