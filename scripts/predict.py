# scripts/predict.py

import pandas as pd
import joblib

from scripts.preprocess import load_and_process_data

def predict_and_save(x_test_path: str, model_path: str, output_path: str):
    X_test_raw = pd.read_csv(x_test_path)
    X_test = load_and_process_data(x_test_path)

    print("Chargement du modèle...")
    model = joblib.load(model_path)
    
    print("Prédictions en cours...")
    y_pred = model.predict(X_test)

    submission = pd.DataFrame(y_pred, columns=["Y1", "Y2", "Y3"])
    submission.insert(0, "ID", X_test_raw["ID"])
    
    submission.to_csv(output_path, index=False)
    print(f"Soumission enregistrée dans {output_path}")

if __name__ == '__main__':
    predict_and_save('data/raw/X_test.csv', 'models/admet_model.pkl', 'submissions/submission.csv')
