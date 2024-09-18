import os
import csv
import json
import torch
import pickle
import optuna
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from scipy.stats import sem
from optuna.samplers import TPESampler
from transformers import AutoModel, AutoTokenizer
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import matthews_corrcoef, make_scorer, roc_auc_score, f1_score
from sklearn.model_selection import ParameterGrid, PredefinedSplit
from optuna.distributions import IntDistribution, CategoricalDistribution
from optuna.integration.sklearn import OptunaSearchCV
from optuna.visualization import plot_param_importances


tokenizer_model = "./model/esm1b_t33_650M_UR50S"
base_model = "./model/non_viral/checkpoint-22836"

string_organism = "Mononegavirales"
organism = "Mononegavirales"

path_results = "./results/"
path_models = "./models/"

method = "nptranfer"
method_name = "NPTranfer"

# Generating fasta path for all folds
train_fasta_folds = ["fold2", "fold3", "fold4", "fold5"]
test_fasta_fold = "fold1"
train_fasta_files = [f"./fasta/folds/{organism}/{fold}.fasta" for fold in train_fasta_folds]
test_fasta_file = f"./fasta/folds/{organism}/{test_fasta_fold}.fasta"

# Fold approach
fold_base_path = f"./input/organismos/taxons/{organism}/"
train_fold_paths = [f"{fold_base_path}fold{i}.csv" for i in [2, 3, 4, 5]]

test_fold = "fold1"
test_path = f"{fold_base_path}{test_fold}.csv"

path_predictions = f"./predictions/{method}/{organism}/{test_fold}.csv"
fold_numbers = 5  # Number of train folds + test_fold

max_seq_length = 512
overlap = 256

def predict_fasta_sequences(fasta_file):
    model = AutoModel.from_pretrained(base_model)
    tokenizer = AutoTokenizer.from_pretrained(tokenizer_model, do_lower_case=False)
    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    sequence_dict = {}

    model.eval()  # Disable dropout

    for record in fasta_sequences:
        sequence = str(record.seq)
        sequence_length = len(sequence)
        embeddings = None

        if sequence_length > max_seq_length:
            # Split the sequence into chunks with overlap
            for i in range(0, sequence_length, max_seq_length - overlap):
                chunk = sequence[i:i + max_seq_length]
                if len(chunk) < max_seq_length:
                    chunk = chunk + ' ' * (max_seq_length - len(chunk))  # Padding if less than max_seq_length
                features = tokenizer.batch_encode_plus(
                    [chunk],
                    add_special_tokens=True,
                    padding='max_length',
                    max_length=max_seq_length,
                    truncation=True,
                    return_tensors='pt',
                    return_attention_mask=True
                )

                with torch.no_grad():
                    outputs = model(features['input_ids'], features['attention_mask'])

                last_hidden_state = outputs[0]
                array_last_hidden_state = last_hidden_state.detach().numpy()

                if embeddings is None:
                    embeddings = array_last_hidden_state[0]
                else:
                    embeddings = np.concatenate((embeddings, array_last_hidden_state[0][overlap:]))
        else:
            # Process the entire sequence at once
            features = tokenizer.batch_encode_plus(
                [sequence],
                add_special_tokens=True,
                padding='max_length',
                max_length=max_seq_length,
                truncation=True,
                return_tensors='pt',
                return_attention_mask=True
            )

            with torch.no_grad():
                outputs = model(features['input_ids'], features['attention_mask'])

            last_hidden_state = outputs[0]
            embeddings = last_hidden_state.detach().numpy()[0]

        sequence_dict[record.id] = {
            'completed_sequence': sequence,
            'embeddings': embeddings
        }

    return sequence_dict

def select_labeled_regions(dataset, sequence_dict):
    list_amino_features = []
    list_amino_label = []

    for idx, row in dataset.iterrows():
        protein_id = row["Info_protein_id"]

        if protein_id not in sequence_dict:
            print(f"{protein_id} is not in sequence_dict")
            continue

        data = sequence_dict[protein_id]
        sequence_slice = data['embeddings']
        start = row["Info_start_pos"]
        end = row["Info_end_pos"]
        label = row["Class"]
        sequence = data['completed_sequence']

        for index_amino in range(start - 1, end):
            list_amino_features.append(sequence_slice[index_amino])
            list_amino_label.append(label)

    return list_amino_features, list_amino_label

# Function to combine the results of multiple calls to predict_fasta_sequences
def predict_fasta_sequences_multiple(files):
    combined_sequence_dict = {}
    for fasta_file in files:
        print(fasta_file)
        sequence_dict = predict_fasta_sequences(fasta_file)
        combined_sequence_dict.update(sequence_dict)
    return combined_sequence_dict

# Call the predict_fasta_sequences function for each fasta file and combine the results
train_sequence_dict = predict_fasta_sequences_multiple(train_fasta_files)
test_sequence_dict = predict_fasta_sequences(test_fasta_file)

# Process each training fold
for i, fold_path in enumerate(train_fold_paths, start=1):
    df_fold = pd.read_csv(fold_path)
    processed_data = select_labeled_regions(df_fold, train_sequence_dict)

    # Convert lists to DataFrames
    features_df = pd.DataFrame(processed_data[0])
    labels_df = pd.DataFrame(processed_data[1])

    # Save features and labels of each fold to CSV files
    features_csv_path = f"{fold_base_path}{method}/processed_fold{i}_features.csv"
    labels_csv_path = f"{fold_base_path}{method}/processed_fold{i}_labels.csv"

    os.makedirs(os.path.dirname(features_csv_path), exist_ok=True)
    os.makedirs(os.path.dirname(labels_csv_path), exist_ok=True)

    features_df.to_csv(features_csv_path, index=False)
    labels_df.to_csv(labels_csv_path, index=False)

    print(f"Processed fold{i} features saved to {features_csv_path}")
    print(f"Processed fold{i} labels saved to {labels_csv_path}")

# Process the test dataset
df_test = pd.read_csv(test_path)
dados_processados = select_labeled_regions(df_test, test_sequence_dict)

# Save features and labels to CSV files
features_csv_path = f"{fold_base_path}{method}/processed_test_features.csv"
labels_csv_path = f"{fold_base_path}{method}/processed_test_labels.csv"
pd.DataFrame(dados_processados[0]).to_csv(features_csv_path, index=False)
pd.DataFrame({"Label": dados_processados[1]}).to_csv(labels_csv_path, index=False)

print(f"Processed test features saved to {features_csv_path}")
print(f"Processed test labels saved to {labels_csv_path}")

loaded_data = {}

# Load the datasets
for i in range(1, fold_numbers):
    features_csv_path = f"{fold_base_path}{method}/processed_fold{i}_features.csv"
    labels_csv_path = f"{fold_base_path}{method}/processed_fold{i}_labels.csv"

    loaded_data[f"fold{i}"] = {
        "features": pd.read_csv(features_csv_path),
        "labels": pd.read_csv(labels_csv_path)
    }

# Load specific test data
test_features_csv_path = f"{fold_base_path}{method}/processed_test_features.csv"
test_labels_csv_path = f"{fold_base_path}{method}/processed_test_labels.csv"

loaded_data["test"] = {
    "features": pd.read_csv(test_features_csv_path),
    "labels": pd.read_csv(test_labels_csv_path)
}


def objective(trial, train_features, train_labels, test_features, test_labels):
    params = {
        'n_estimators': trial.suggest_int('n_estimators', 100, 500),
        'max_depth': trial.suggest_categorical('max_depth', [None, 10, 20, 30]),
        'min_samples_split': trial.suggest_int('min_samples_split', 2, 10),
        'min_samples_leaf': trial.suggest_int('min_samples_leaf', 1, 10),
        'max_features': trial.suggest_categorical('max_features', ['sqrt', 'log2', None]),
        'bootstrap': trial.suggest_categorical('bootstrap', [True, False]),
        'criterion': trial.suggest_categorical('criterion', ['gini', 'entropy']),
        'random_state': 42
    }

    model = RandomForestClassifier(**params, n_jobs=-1)
    model.fit(train_features, train_labels)

    train_probs = model.predict_proba(train_features)[:, 1]
    thresholds = np.arange(0, 1, 0.001)
    mcc_scores = [matthews_corrcoef(train_labels, (train_probs > t).astype(int)) for t in thresholds]
    max_mcc_index = np.argmax(mcc_scores)
    threshold = thresholds[max_mcc_index]

    test_probs = model.predict_proba(test_features)[:, 1]
    predicted_labels = (test_probs > threshold).astype(int)
    test_mcc = matthews_corrcoef(test_labels, predicted_labels)

    # Save model, threshold and params in trial for later use
    trial.set_user_attr('model', model)
    trial.set_user_attr('threshold', threshold)
    trial.set_user_attr('params', params)

    return test_mcc

def perform_validation_holdout(loaded_data):
    train_features = loaded_data["fold1"]["features"]
    train_labels = loaded_data["fold1"]["labels"].iloc[:, 0]
    test_features = loaded_data["test"]["features"]
    test_labels = loaded_data["test"]["labels"].iloc[:, 0]

    study = optuna.create_study(direction='maximize')
    study.optimize(lambda trial: objective(trial, train_features, train_labels, test_features, test_labels), n_trials=200)

    best_trial = study.best_trial
    final_model = best_trial.user_attrs['model']
    mcc_best_threshold = best_trial.user_attrs['threshold']
    best_hyperparams = best_trial.params

    # Get prediction probabilities for all data
    all_probs = final_model.predict_proba(test_features)[:, 1]

    f1_best_threshold = 0.5

    # Save the best model
    model_path = os.path.join(path_models, f"{method}_{organism}.pkl")
    with open(model_path, 'wb') as f:
        pickle.dump(final_model, f)
    print(f"Best model saved to {model_path}")

    return final_model, mcc_best_threshold, f1_best_threshold, best_hyperparams

# Define custom MCC scoring function
def mcc_scorer(y_true, y_pred):
    return matthews_corrcoef(y_true, y_pred)

mcc_scorer = make_scorer(mcc_scorer)

# Define custom AUC-ROC scoring function
def auc_roc_scorer(y_true, y_pred):
    return roc_auc_score(y_true, y_pred)

roc_auc_scorer = make_scorer(auc_roc_scorer, needs_proba=True)

def perform_cross_validation(loaded_data):
    # Concatenate all folds into a single dataset
    all_features = pd.concat([loaded_data[f"fold{i}"]["features"] for i in range(1, len(loaded_data))], ignore_index=True)
    all_labels = pd.concat([loaded_data[f"fold{i}"]["labels"] for i in range(1, len(loaded_data))], ignore_index=True).iloc[:, 0]

    # Initialize lists to store training and validation indices
    train_indices_list = []
    val_indices_list = []

    # Iterate over the folds in the loaded_data dictionary to define training and validation indices
    for i in range(1, len(loaded_data)):
        train_indices = []
        val_indices = loaded_data[f"fold{i}"]["features"].index.tolist()

        validation_labels = loaded_data[f"fold{i}"]["labels"].iloc[:, 0]
        for j in range(1, len(loaded_data)):
            if j != i:
                train_indices.extend(loaded_data[f"fold{j}"]["features"].index.tolist())
        train_indices_list.append(train_indices)
        val_indices_list.append(val_indices)

    # Define a custom class for StratifiedKFold
    class CustomStratifiedKFold:
        def __init__(self, train_indices_list, val_indices_list):
            self.train_indices_list = train_indices_list
            self.val_indices_list = val_indices_list
            self.n_splits = len(train_indices_list)

        def split(self, X, y=None, groups=None):
            for i in range(self.n_splits):
                yield self.train_indices_list[i], self.val_indices_list[i]

        def get_n_splits(self, X=None, y=None, groups=None):
            return self.n_splits

    custom_skf = CustomStratifiedKFold(train_indices_list, val_indices_list)

    # Hyperparameters to be optimized by OptunaSearchCV
    param_distributions = {
        'n_estimators': IntDistribution(low=100, high=500),
        'max_depth': CategoricalDistribution(choices=[None, 10, 20, 30]),
        'min_samples_split': IntDistribution(low=2, high=10),
        'min_samples_leaf': IntDistribution(low=1, high=10),
        'max_features': CategoricalDistribution(choices=['sqrt', 'log2', None]),
        'bootstrap': CategoricalDistribution(choices=[True, False]),
        'criterion': CategoricalDistribution(choices=['gini', 'entropy'])
    }

    model = RandomForestClassifier(random_state=42, n_jobs=-1)
    optuna_search = OptunaSearchCV(
        estimator=model,
        param_distributions=param_distributions,
        n_trials=100,
        random_state=42,
        cv=custom_skf,
        n_jobs=-1,
        scoring=roc_auc_scorer
    )

    optuna_search.fit(all_features, all_labels)

    # Obtain the best hyperparameters
    best_hyperparams = optuna_search.best_params_
    final_model = RandomForestClassifier(**best_hyperparams, random_state=42, n_jobs=-1)

    final_model.fit(all_features, all_labels)

    all_probs = final_model.predict_proba(all_features)[:, 1]

    # Determine the best threshold maximizing the MCC
    thresholds = np.arange(0, 1, 0.001)
    mcc_scores = [matthews_corrcoef(all_labels, all_probs >= t) for t in thresholds]
    mcc_best_threshold_index = np.argmax(mcc_scores)
    mcc_best_threshold = thresholds[mcc_best_threshold_index]

    f1_scores = [f1_score(all_labels, all_probs >= t) for t in thresholds]
    f1_best_threshold_index = np.argmax(f1_scores)
    f1_best_threshold = thresholds[f1_best_threshold_index]

    # Save the best model
    model_path = os.path.join(path_models, f"{method}_{organism}.pkl")
    with open(model_path, 'wb') as f:
        pickle.dump(final_model, f)
    print(f"Best model saved to {model_path}")

    return final_model, mcc_best_threshold, f1_best_threshold, best_hyperparams

def evaluate_model_on_test(test_features, test_labels, model, mcc_threshold, f1_threshold):
    test_probs = model.predict_proba(test_features)[:, 1]

    # Final predictions for each threshold
    predictions_mcc = test_probs >= mcc_threshold
    predictions_f1 = test_probs >= f1_threshold

    # Calculate metrics based on final predictions
    auc = roc_auc_score(test_labels, test_probs)

    # Metrics for MCC Threshold
    f1_at_mcc_threshold = f1_score(test_labels, predictions_mcc)
    mcc_at_mcc_threshold = matthews_corrcoef(test_labels, predictions_mcc)

    # Metrics for F1 Threshold
    f1_at_f1_threshold = f1_score(test_labels, predictions_f1)
    mcc_at_f1_threshold = matthews_corrcoef(test_labels, predictions_f1)

    print(f"AUC: {auc}")
    print(f"At MCC Threshold -> F1 Score: {f1_at_mcc_threshold}, MCC: {mcc_at_mcc_threshold}")

    # Prepare data for CSV
    results = pd.DataFrame([{
        "Method": method,
        "Organism": organism,
        "Test probabilities": str(test_probs.tolist()),
        "MCC threshold": mcc_threshold,
        "F1 threshold": f1_threshold
    }])

    # Path to results file
    results_file_path = os.path.join(path_results, f"{method}_{organism}_test_prob_thresholds.csv")

    # Check if the file already exists and load the existing content if necessary
    if os.path.exists(results_file_path):
        existing_results = pd.read_csv(results_file_path)
        results = pd.concat([existing_results, results], ignore_index=True)

    # Save the results to CSV
    results.to_csv(results_file_path, index=False)

    return {
        'auc': auc,
        'mcc_threshold': {'f1': f1_at_mcc_threshold, 'mcc': mcc_at_mcc_threshold},
        'f1_threshold': {'f1': f1_at_f1_threshold, 'mcc': mcc_at_f1_threshold}
    }

def main(loaded_data, organism, method, path_results):
    if len(loaded_data) > 2:
        best_model, mcc_best_threshold, f1_best_threshold, best_hyperparams = perform_cross_validation(loaded_data)

        # Load test data
        test_features = loaded_data["test"]["features"]
        test_labels = loaded_data["test"]["labels"].iloc[:, 0]

        # Evaluate the model on the test set
        results_dict = evaluate_model_on_test(test_features, test_labels, best_model, mcc_best_threshold, f1_best_threshold)

        # Adding results for MCC threshold
        results = pd.DataFrame({
            "Dataset": [organism] * 3,
            "Method": [method] * 3,
            "Metric": ["AUC", "F1 (MCC Threshold)", "MCC (MCC Threshold)"],
            "Value": [results_dict['auc'], results_dict['mcc_threshold']['f1'], results_dict['mcc_threshold']['mcc']],
            "Best threshold": ["", mcc_best_threshold, mcc_best_threshold]
        })

        # Concatenating the results
        results["Value"] = results["Value"].round(3)

        # Including the best hyperparameters as a single string
        best_hyperparams_str = ', '.join([f"{key}: {value}" for key, value in best_hyperparams.items()])
        results["Best hyperparameters"] = [""] * len(results)
        results.loc[0, "Best hyperparameters"] = best_hyperparams_str

        # Path to results file
        results_file_path = os.path.join(path_results, "test_metrics_summary.csv")

        # Checking if the file already exists and loading the existing content if necessary
        if os.path.exists(results_file_path):
            existing_results = pd.read_csv(results_file_path)
            results = pd.concat([existing_results, results], ignore_index=True)

        # Saving the results to CSV
        results.to_csv(results_file_path, index=False)

        print(f"Results saved to {results_file_path}")
    else:
        best_model, mcc_best_threshold, f1_best_threshold, best_hyperparams = perform_validation_holdout(loaded_data)

        # Load test data
        test_features = loaded_data["test"]["features"]
        test_labels = loaded_data["test"]["labels"].iloc[:, 0]

        # Evaluate the model on the test set
        results_dict = evaluate_model_on_test(test_features, test_labels, best_model, mcc_best_threshold, f1_best_threshold)

        # Adding results for MCC threshold
        results = pd.DataFrame({
            "Dataset": [organism] * 3,
            "Method": [method] * 3,
            "Metric": ["AUC", "F1 (MCC Threshold)", "MCC (MCC Threshold)"],
            "Value": [results_dict['auc'], results_dict['mcc_threshold']['f1'], results_dict['mcc_threshold']['mcc']],
            "Best threshold": ["", mcc_best_threshold, mcc_best_threshold]
        })

        # Adding results for F1 threshold
        additional_results = pd.DataFrame({
            "Dataset": [organism] * 2,
            "Method": [method] * 2,
            "Metric": ["F1 (F1 Threshold)", "MCC (F1 Threshold)"],
            "Value": [results_dict['f1_threshold']['f1'], results_dict['f1_threshold']['mcc']],
            "Best threshold": [f1_best_threshold, f1_best_threshold]
        })

        # Concatenating the results
        results = pd.concat([results, additional_results], ignore_index=True)

        results["Value"] = results["Value"].round(3)

        # Including the best hyperparameters as a single string
        best_hyperparams_str = ', '.join([f"{key}: {value}" for key, value in best_hyperparams.items()])
        results["Best hyperparameters"] = [""] * len(results)
        results.loc[0, "Best hyperparameters"] = best_hyperparams_str

        # Path to results file
        results_file_path = os.path.join(path_results, "test_metrics_summary.csv")

        # Checking if the file already exists and loading the existing content if necessary
        if os.path.exists(results_file_path):
            existing_results = pd.read_csv(results_file_path)
            results = pd.concat([existing_results, results], ignore_index=True)

        # Saving the results to CSV
        results.to_csv(results_file_path, index=False)

        print(f"Results saved to {results_file_path}")

main(loaded_data, organism, method, path_results)

# Path of the saved model
model_path = os.path.join(path_models, f"{method}_{organism}.pkl")

# Load the saved model
with open(model_path, 'rb') as f:
    model = pickle.load(f)
print(f"Model loaded from {model_path}")

# Function to generate the CSV
def generate_csv(sequence_dict, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Protein_id', 'amino', 'Predicted_prob', 'Position'])

        for protein_id, data in sequence_dict.items():
            sequence = data['completed_sequence']
            predictions = data['embeddings']
            sequence_length = len(sequence)
            for i in range(sequence_length):
                amino_acid = sequence[i]  # Get the amino acid from the sequence
                predicted_prob = predictions[i]  # Get the prediction for the specific position
                writer.writerow([protein_id, amino_acid, predicted_prob, i + 1])

    print(f'CSV file {output_file} generated successfully.')

# Load the saved model
with open(model_path, 'rb') as f:
    model = pickle.load(f)
print(f"Model loaded from {model_path}")

# Generate predictions for the complete .fasta using the loaded model
for protein_id, data in test_sequence_dict.items():
    features = np.array(data['embeddings'])
    predictions = model.predict_proba(features)[:, 1]  # RandomForest returns probabilities for each class
    data['embeddings'] = predictions

# Define the output file path for the CSV
output_file = f"./predictions/{method}/{organism}/{test_fold}/formatted_raw_output.csv"
generate_csv(test_sequence_dict, output_file)

# Generate the CSV for the test_fold
file_path = f'{path_results}test_metrics_summary.csv'
df = pd.read_csv(file_path)

# Sorting the DataFrame
df_sorted = df.sort_values(by=['Metric', 'Dataset'])

# Save the sorted DataFrame to a new CSV file
sorted_file_path = f'{path_results}sorted_metrics_summary.csv'
df_sorted.to_csv(sorted_file_path, index=False)

# Identifying the best method for each dataset by Metric
best_methods = df_sorted.loc[df_sorted.groupby(['Dataset', 'Metric'])['Value'].idxmax()]

best_methods = best_methods.drop(columns=['Best threshold', 'Best hyperparameters'])

best_methods.rename(columns={'Method': 'Best Method'}, inplace=True)

# Save the best methods to a single CSV file
best_methods_file_path = f'{path_results}best_methods_summary.csv'
best_methods.to_csv(best_methods_file_path, index=False)

df = pd.read_csv(best_methods_file_path)

# Count the number of times each method appears as the best method
method_counts = df['Best Method'].value_counts()

# Determine the best overall method (the method with the highest count)
best_overall_method = method_counts.idxmax()

# Print the counts and the best overall method
print("Method counts:\n", method_counts)
print("\nBest overall method:", best_overall_method)
