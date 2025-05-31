import sys
import pandas as pd
from joblib import load
from sklearn.metrics import matthews_corrcoef, f1_score, roc_auc_score, balanced_accuracy_score, precision_score, recall_score, confusion_matrix
import os

if __name__ == "__main__":
    taxa = ["bpertussis", "corynebacterium", "orthopoxvirus", "ecoli", "enterobacteriaceae", "lentivirus", "mtuberculosis", "paeruginosa", "smansoni", "tgondii", "pfalciparum", "ctrachomatis", "human_gammaherpesvirus_4", "influenza_a", "cdifficile", "filoviridae", "measles_morbilivirus", "ovolvulus", "mononegavirales"]
    taxa_string = ["B. pertussis", "Corynebacterium", "Orthopoxvirus", "E. coli", "Enterobacteriaceae", "Lentivirus", "M. tuberculosis", "P. aeruginosa", "S. mansoni", "T. gondii", "P. falciparum", "C. trachomatis", "Human Gammaherpesvirus 4", "Influenza A", "C. difficile", "Filoviridae", "Measles morbilivirus", "Ovolvulus", "Mononegavirales"]
    thresholds = [0.422, 0.5, 0.5, 0.492, 0.432, 0.6, 0.507, 0.445, 0.262, 0.443, 0.571, 0.516, 0.435, 0.559, 0.241, 0.365, 0.5, 0.501, 0.325]

    print("Diretório de trabalho:", os.getcwd())
    
    taxa_inputs = sys.argv[1:]  
    if "all" in taxa_inputs:
        taxa_to_process = taxa
    else:
        taxa_to_process = [t for t in taxa_inputs if t in taxa]

    for taxa_input in taxa_to_process:
        index = taxa.index(taxa_input)
        string_taxa = taxa_string[index]
        best_threshold = thresholds[index]

        
        model_path = f"./models/epitopetransfer_{taxa_input}.pkl"
        model = load(model_path)

        test_features_csv_path = f"./input/epitopetransfer/{taxa_input}/processed_test_features.csv"
        test_labels_csv_path = f"./input/epitopetransfer/{taxa_input}/test_labels.csv"

        test_features = pd.read_csv(test_features_csv_path)
        test_labels = pd.read_csv(test_labels_csv_path).iloc[:, 0]

        test_probs = model.predict_proba(test_features)[:, 1]
        test_predictions = (test_probs > best_threshold).astype(int)

        # metrics
        test_mcc = matthews_corrcoef(test_labels, test_predictions)
        f1 = f1_score(test_labels, test_predictions)
        auc_score = roc_auc_score(test_labels, test_probs)
        bacc = balanced_accuracy_score(test_labels, test_predictions)
        ppv = precision_score(test_labels, test_predictions, zero_division=0)
        sens = recall_score(test_labels, test_predictions, zero_division=0) 

        tn, fp, fn, tp = confusion_matrix(test_labels, test_predictions).ravel()
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        npv = tn / (tn + fn) if (tn + fn) > 0 else 0

        title = f"{string_taxa}"
        results = (
            f"AUC: {auc_score:.3f}\n"
            f"F1 : {f1:.3f}\n"
            f"MCC: {test_mcc:.3f}\n"
            f"BACC: {bacc:.3f}\n"
            f"PPV: {ppv:.3f}\n"
            f"NPV: {npv:.3f}\n"
            f"SENS: {sens:.3f}\n"
            f"SPEC: {specificity:.3f}"
        )
        print(f"{title}\n{'-' * len(title)}\n{results}\n")

