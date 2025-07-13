# scripts/utils.py

import os

def ensure_dir(path: str):
    """Crée un dossier s'il n'existe pas."""
    if not os.path.exists(path):
        os.makedirs(path)
