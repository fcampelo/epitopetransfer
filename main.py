import sys
import pandas as pd
from joblib import load
from sklearn.metrics import matthews_corrcoef, f1_score, roc_auc_score

if __name__ == "__main__":
    taxa = ["bpertussis", "corynebacterium"]
    taxa_string = ["B. pertussis", "Corynebacterium"]
    thresholds = [0.801, 0.468]
    method = "epitopetransfer"
    
    taxa_inputs = sys.argv[1:]  # Capture all arguments 
    if "all" in taxa_inputs:
        taxa_to_process = taxa
    else:
        taxa_to_process = [t for t in taxa_inputs if t in taxa]

    for taxa_input in taxa_to_process:
        index = taxa.index(taxa_input)
        string_taxa = taxa_string[index]
        best_threshold = thresholds[index]

        fold_base_path = "./results/"
        

        model_path = f"./models/{taxa_input}.joblib"
        model = load(model_path)

        test_features_csv_path = f"./input/{method}/processed_test_features.csv"
        test_labels_csv_path = f"./input/{method}/processed_test_labels.csv"

        test_features = pd.read_csv(test_features_csv_path)
        test_labels = pd.read_csv(test_labels_csv_path).iloc[:, 0]

        test_probs = model.predict_proba(test_features)[:, 1]
        test_predictions = (test_probs > best_threshold).astype(int)

        test_mcc = matthews_corrcoef(test_labels, test_predictions)
        f1 = f1_score(test_labels, test_predictions, average='macro')
        auc_score = roc_auc_score(test_labels, test_probs)

        # Print the evaluation results
        title = f"{string_taxa}"
        results = f"AUC: {auc_score:.3f}\nF1 : {f1:.3f}\nMCC: {test_mcc:.3f}"
        print(f"{title}\n{'-' * len(title)}\n{results}\n")

