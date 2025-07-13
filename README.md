# 🧬 Prédiction des Profils ADMET par Apprentissage Automatique

## 📌 Contexte

Un médicament efficace ne se limite pas seulement à sa capacité d'interagir avec une cible biologique spécifique ; ses propriétés ADMET (Absorption, Distribution, Métabolisme, Excrétion, Toxicité) sont cruciales pour déterminer sa sécurité et son efficacité réelle chez les patients. Cependant, l'évaluation expérimentale de ces propriétés est coûteuse et chronophage.

L'apprentissage automatique offre une solution prometteuse en permettant de prédire ces propriétés rapidement et efficacement, accélérant ainsi la découverte de nouveaux médicaments tout en réduisant les coûts et les risques associés.

## 🎯 Objectif du Projet

Ce projet vise à développer un modèle de classification multi-étiquettes supervisée capable de prédire trois propriétés ADMET pour diverses molécules. Chaque propriété est indiquée par une étiquette binaire (1 = plage optimale, 0 = hors plage).

Les prédictions réalisées permettront d'identifier rapidement des candidats-médicaments possédant des profils ADMET favorables.

## 📂 Description des Données

Vous disposez de cinq fichiers principaux :

* `X_train.csv` : 1940 molécules avec leurs identifiants uniques et représentations SMILES.
* `y_train.csv` : Étiquettes ADMET pour chaque molécule d'entraînement.
* `X_test.csv` : 486 molécules à prédire.
* `random_submission_example.csv` : Exemple du format de soumission attendu.
* `supplementary_files` : Notebook Jupyter introductif présentant le défi et un modèle de référence.

### Structure des données :

| Colonne  | Description                                   |
| -------- | --------------------------------------------- |
| `ID`     | Identifiant unique des molécules              |
| `SMILES` | Chaîne de caractères représentant la molécule |
| `Y1`     | Étiquette pour la propriété ADMET 1           |
| `Y2`     | Étiquette pour la propriété ADMET 2           |
| `Y3`     | Étiquette pour la propriété ADMET 3           |


### Structure du projet 

Projet_ADMET/
├── data/
│   ├── raw/                 # Données brutes (X_train.csv, X_test.csv, y_train.csv)
│   └── processed/           # Données nettoyées ou vectorisées
│
├── notebooks/
│   ├── exploration.ipynb        # Analyse exploratoire initiale
│   ├── preprocessing.ipynb      # Transformation des SMILES, nettoyage, empreintes
│   ├── baseline_model.ipynb     # Modèle de référence
│   ├── modeling.ipynb           # Optimisation des modèles
│   └── evaluation.ipynb         # Évaluation et analyse des résultats
│
├── scripts/
│   ├── preprocess.py            # Prétraitement et vectorisation
│   ├── train.py                 # Entraînement des modèles
│   ├── predict.py               # Prédictions sur les données de test
│   └── utils.py                 # Fonctions utilitaires communes
│
├── app/
│   ├── streamlit_app.py         # Application Streamlit principale
│   ├── components/              # Composants personnalisés (visualisation moléculaire, sliders, etc.)
│   └── assets/                  # Images, logos, fichiers statiques
│
├── models/                      # Modèles sauvegardés (.pkl, .joblib)
│
├── submissions/                 # Fichiers de prédiction au format d’envoi
│
├── reports/                     # Graphiques, métriques, logs, courbes d'apprentissage
│
├── requirements.txt             # Dépendances (inclure streamlit, rdkit, scikit-learn, etc.)
├── README.md                    # Présentation du projet
└── .gitignore                   # Fichiers à exclure du contrôle de version

## 🛠️ Méthodologie

1. **Prétraitement des données** :

   * Conversion des représentations SMILES en empreintes moléculaires (vecteurs binaires).

2. **Modélisation** :

   * Application de techniques d'apprentissage automatique (arbres de décision, forêts aléatoires, modèles ensemblistes, réseaux de neurones).

3. **Évaluation et Optimisation** :

   * Métrique utilisée : **Score F1 micro-moyenné** (combinant précision et rappel globalement).

## 📈 Métrique d'Évaluation

Le modèle sera évalué avec le **score F1 micro-moyenné**, défini comme suit :

$$
Precision_{micro} = \frac{\sum_c TP_c}{\sum_c (TP_c + FP_c)}
$$

$$
Recall_{micro} = \frac{\sum_c TP_c}{\sum_c (TP_c + FN_c)}
$$

$$
F1_{micro} = 2 \times \frac{Precision_{micro} \times Recall_{micro}}{Precision_{micro} + Recall_{micro}}
$$


## 🚀 Soumission des Résultats

Vos prédictions doivent être soumises au format du fichier `random_submission_example.csv`.

## 🧑‍💻 Contributeurs

Oumar Abdramane ALLAWAN
