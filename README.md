# EpitopeTransfer

This repository contains data and code for the paper EpitopeTransfer: Improving linear B-cell epitope prediction through phylogeny-aware transfer learning by Lindeberg Pessoa Leite, Teófilo E. de Campos, Francisco Pereira Lobo and Felipe Campelo.
s
## Table of Contents
- [Dependencies](#dependencies)
- [Run Models for Paper Metrics](#run-models-for-paper-metrics)
- [Make full-pipeline predictions](#make-full-pipeline-predictions)
- [Citation](#citation)
- [Contact](#contact)

## Dependencies

This project requires Python 3.10. Follow the instructions below to set up your environment and install the required dependencies.


1. **Install Python 3.10**

   Ensure that Python 3.10 is installed on your system. If not, you can download it from the [official Python website](https://www.python.org/downloads/release/python-3100/) or install it using your system's package manager. For Ubuntu-based distributions, you can use:

   ```bash
   sudo apt update
   sudo apt install python3.10

For those who prefer Conda, create a project environment following the [Conda environment management guide](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

2. ** Project dependencies**

  ```bash
  pip install -r requirements.txt


## Run Models for Paper Metrics

To run the models for specific taxa or for all taxa included in the study, use the following command format in the terminal:

```bash
python main.py [taxa]

Replace [taxa] with the name of the taxa you wish to process from the list below, or use all to process all available taxa. Example:

```bash
python main.py bpertussis

### Available Taxa

- `bpertussis` - For running the model on B. pertussis data.
- `corynebacterium` - For running the model on Corynebacterium data.
- `all` - Use this option to run the models for all the above taxa.


## Make full-pipeline predictions

Run EpitopeTransfer model with your own data

## Citation
If you find our work useful in your research, please consider citing:
```bibtex
@inproceedings{lindeberg2024EpitopeTransfer,
  title={EpitopeTransfer: Improving linear B-cell epitope prediction through phylogeny-aware transfer learning},
  author={Lindeberg Leite and Teófilo de Campos and Francisco Lobo and Felipe Campelo},
  booktitle={Briefings in Bioinformatics},
  year={2024},
}
````
## Contact
Email: [f.campelo@aston.ac.uk](mailto:f.campelo@aston.ac.uk) - Principal investigator <br>
Email: [lindpessoa@gmail.com](mailto:lindpessoa@gmail.com)
