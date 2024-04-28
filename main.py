import sys
import pandas as pd
from joblib import load
from sklearn.metrics import matthews_corrcoef, f1_score, roc_auc_score

if __name__ == "__main__":
    taxa = ["bpertussis", "corynebacterium", "orthopoxvirus", "ecoli", "enterobacteriaceae", "lentivirus", "mtuberculosis", "paeruginosa", "sars_cov2", "smansoni", "tgondii", "pfalciparum"]
    taxa_string = ["B. pertussis", "Corynebacterium", "Orthopoxvirus", "E. coli", "Enterobacteriaceae", "Lentivirus", "M. tuberculosis", "P. aeruginosa", "SARS-Cov-2", "S. mansoni", "T. gondii", "P. falciparum"]
    thresholds = [0.801, 0.468, 0.497, 0.481, 0.201, 0.841, 0.226, 0.7, 0.169, 0.274, 0.69, 0.879]

    
    taxa_inputs = sys.argv[1:]  
    if "all" in taxa_inputs:
        taxa_to_process = taxa
    else:
        taxa_to_process = [t for t in taxa_inputs if t in taxa]

    for taxa_input in taxa_to_process:
        index = taxa.index(taxa_input)
        string_taxa = taxa_string[index]
        best_threshold = thresholds[index]

        
        model_path = f"./models/{taxa_input}.pkl"
        model = load(model_path)

        test_features_csv_path = f"./input/{taxa_input}/processed_test_features.csv"
        test_labels_csv_path = f"./input/{taxa_input}/test_labels.csv"

        test_features = pd.read_csv(test_features_csv_path)
        test_labels = pd.read_csv(test_labels_csv_path).iloc[:, 0]

        test_probs = model.predict_proba(test_features)[:, 1]
        test_predictions = (test_probs > best_threshold).astype(int)

        test_mcc = matthews_corrcoef(test_labels, test_predictions)
        f1 = f1_score(test_labels, test_predictions)
        auc_score = roc_auc_score(test_labels, test_probs)

        # Print the evaluation results
        title = f"{string_taxa}"
        results = f"AUC: {auc_score:.3f}\nF1 : {f1:.3f}\nMCC: {test_mcc:.3f}"
        print(f"{title}\n{'-' * len(title)}\n{results}\n")

