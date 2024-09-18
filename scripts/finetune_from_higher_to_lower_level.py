
import logging
from transformers import AutoModelForSequenceClassification
from transformers import AutoTokenizer
import torch
import pandas as pd
from transformers import BertTokenizer, BertForSequenceClassification
from transformers import TrainingArguments, Trainer
from sklearn.metrics import accuracy_score, recall_score, precision_score, f1_score
from transformers import AutoModelForSequenceClassification
import numpy as np
from sklearn.metrics import accuracy_score, recall_score, precision_score, f1_score, matthews_corrcoef


#  This code performs fine-tuning of a model from higher taxonomic levels 
#  to lower taxonomic levels


logging.basicConfig(level=logging.DEBUG)

# Replace here with your paths
tokenizer_model = "./model/esm1b_t33_650M_UR50S"
base_model = "./model/esm1b_t33_650M_UR50S"

input_training_data = "./input/organismos/taxons/Mononegavirales/train_highest_minus_lowest.csv"
input_test_data = "./input/organismos/taxons/Mononegavirales/test_highest_minus_lowest.csv"
generated_model = "./model/mononegavirales_highest_minus_lowest"


tokenizer = AutoTokenizer.from_pretrained(tokenizer_model, do_lower_case=False)
model = AutoModelForSequenceClassification.from_pretrained(base_model, num_labels=2)



# Create torch dataset
class Dataset(torch.utils.data.Dataset):
    def __init__(self, encodings, labels=None):
        self.encodings = encodings
        self.labels = labels

    def __getitem__(self, idx):
        item = {key: torch.tensor(val[idx]) for key, val in self.encodings.items()}
        if self.labels:
            item["labels"] = torch.tensor(self.labels[idx])
        return item

    def __len__(self):
        return len(self.encodings["input_ids"])


# Passing throuth the label region considering the window size on the whole protein
def segmentBycontext(sequences, start, end, classes,  window):

    center = []
    left = []
    right = []
    labels = []
    count=0;
    window_size = window

    for sequence, start, end, label in zip(sequences, start, end, classes):
        length_sequence = len(sequence)
        position=start-1;

        while position < end:

            center_amino = sequence[position]

            if position > window_size:
                left_aminos = sequence[position - window_size:position]

            else:
                left_aminos = sequence[0:position]

            if window_size < (length_sequence - position) :
                right_aminos = sequence[position + 1:position + 1 + window_size]

            else:
                right_aminos = sequence[position + 1:length_sequence]


            center.append(center_amino)
            left.append(left_aminos)
            right.append(right_aminos)
            labels.append(label)



            position = position + 1;
            count = count + 1;


    return center, left, right, labels

def process_data_for_bert(data, tokenizer):

    sentence_center = list(data["center"])
    sentence_left = list(data["left"])
    sentence_right = list(data["right"])

    left_center_right = []

    for (left, center, right) in zip(sentence_left, sentence_center, sentence_right):
        left_center_right.append(left + center + right)

    y = list(data["label"])

    X_tokenized = tokenizer(left_center_right, padding=True, truncation=True, max_length=512)
    dataset = Dataset(X_tokenized, y)

    return dataset

def data_processing_pipeline(data, window_size, tokenizer):

    start = list(data["Info_start_pos"])
    end = list(data["Info_end_pos"])

    sequence = data['Info_sequence'].tolist()
    labels = data['Class'].tolist()

    center, left, right, labels =  segmentBycontext(sequence, start, end, labels, window_size)

    sentences = pd.DataFrame({'center':center, 'left':left, 'right':right, 'label':labels})

    dataset = process_data_for_bert(sentences, tokenizer)


    return dataset


##  Fine-tune pretrained model

batch_size=8

tokenizer = AutoTokenizer.from_pretrained(tokenizer_model, do_lower_case=False)
model = AutoModelForSequenceClassification.from_pretrained(base_model, num_labels=2)

model.cuda()

train_dataset = pd.read_csv(input_training_data)
val_dataset = pd.read_csv(input_test_data)

train_dataset = data_processing_pipeline(train_dataset, 512, tokenizer)
val_dataset = data_processing_pipeline(val_dataset, 512, tokenizer)

# Define Trainer parameters
def compute_metrics(p):
    pred, labels = p
    pred = np.argmax(pred, axis=1)

    accuracy = accuracy_score(y_true=labels, y_pred=pred)
    recall = recall_score(y_true=labels, y_pred=pred)
    precision = precision_score(y_true=labels, y_pred=pred)
    f1 = f1_score(y_true=labels, y_pred=pred)
    mcc = matthews_corrcoef(y_true=labels, y_pred=pred)

    return {"accuracy": accuracy, "precision": precision, "recall": recall, "f1": f1, "mcc": mcc}


args = TrainingArguments(
    generated_model,
    evaluation_strategy = "epoch",
    save_strategy = "epoch",
    learning_rate=1e-5,
    per_device_train_batch_size=batch_size,
    per_device_eval_batch_size=batch_size,
    num_train_epochs=3,
    weight_decay=0.01,
    load_best_model_at_end=True,
    metric_for_best_model='f1',
    optim="adamw_torch"
)

def model_init():
    return model

trainer = Trainer(
    model_init=model_init,
    args=args,
    train_dataset=train_dataset,
    eval_dataset=val_dataset,
    tokenizer=tokenizer,
    compute_metrics=compute_metrics
)

# New code - wrap collator in a dictionary
data_collator = trainer.data_collator
trainer.data_collator = lambda data: dict(data_collator(data))
# End new code

trainer.train()

trainer.save_model(generated_model)

