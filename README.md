# EpitopeTransfer
[![DOI](https://zenodo.org/badge/859842116.svg)](https://zenodo.org/doi/10.5281/zenodo.13811925)

This repository contains data and code for the paper **EpitopeTransfer:
Improving linear B-cell epitope prediction through phylogeny-aware
transfer learning** by Lindeberg Pessoa Leite, Teófilo E. de Campos,
Francisco Pereira Lobo and Felipe Campelo.

## Table of Contents

-   [Dependencies](#dependencies)
-   [Running Models for Published
    Metrics](#running-models-for-paper-metrics)
-   [Running the analyses](#Running-the-analyses)
-   [Contact](#contact)

## Dependencies {#dependencies}

This project requires Python 3.10, although other versions are likely
compatible. Follow the instructions below to set up your environment and
install the required dependencies.

1.  **Install Python 3.10**

Ensure that Python 3.10 is installed on your system. If not, you can download it from the [official Python website](https://www.python.org/downloads/release/python-3100/) or install it using your system's package manager. For Ubuntu-based distributions, you can use:

``` bash
  sudo apt update
  sudo apt install python3.10
```

For those who prefer Conda, create a project environment following the [Conda environment management guide](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

2.  **Project dependencies**

``` bash
  pip install -r requirements.txt
```

## Run Models for Paper Metrics {#run-models-for-paper-metrics}

To run the models for specific taxa or for all taxa included in the
study, use the following command format in the terminal:

``` bash
  python main.py [taxa]
```

Replace [taxa] with the name of the taxa you wish to process from the list below, or use all to process all available taxa. Example:

``` bash
  python main.py bpertussis
```

Available Taxa

| **Taxa**                                        | **Taxa**                                                 |
|-----------------------------------|-------------------------------------|
| **bpertussis**: *Bordetella pertussis*          | **filoviridae**: *Filoviridae*                           |
| **corynebacterium**: *Corynebacterium*          | **ovolvulus**: *Onchocerca volvulus*                     |
| **orthopoxvirus**: *Orthopoxvirus*              | **ctrachomatis**: *Chlamydia trachomatis*                |
| **ecoli**: *Escherichia coli*                   | **human_gammaherpesvirus_4**: *Human Gammaherpesvirus 4* |
| **enterobacteriaceae**: *Enterobacteriaceae*    | **influenza_a**: *Influenza A*                           |
| **lentivirus**: *Lentivirus*                    | **cdifficile**: *Clostridioides difficile*               |
| **mtuberculosis**: *Mycobacterium tuberculosis* | **measles_morbilivirus**: *Measles morbillivirus*        |
| **paeruginosa**: *Pseudomonas aeruginosa*       | **mononegavirales**: *Mononegavirales*                   |
| **sars_cov2**: *SARS-CoV-2*                     |                                                          |
| **smansoni**: *Schistosoma mansoni*             |                                                          |
| **tgondii**: *Toxoplasma gondii*                |                                                          |
| **pfalciparum**: *Plasmodium falciparum*        |                                                          |

obs: To run the taxa in the second column, upgrade scikit-learn to
version 1.3.2 (pip install scikit-learn-1.3.2)

## Running the analyses

The analysis of the results was done using R version 4.4.1 (reproducible
using the script under folder `./R`). The main packages used in the
analysis were:

-   dplyr_1.1.4
-   tidyr_1.3.1
-   yardstick_1.3.1
-   pROC_1.18.5
-   multcomp_1.4-26
-   ggplot2_3.5.1
-   ggrepel_0.9.5
-   see_0.9.0
-   stringr_1.5.1
-   wrappedtools_0.9.5

The full details of the R session used in the analysis are available in `/R/SessionInfo.txt`

*****
### Contact Email:  
[Felipe Campelo](mailto:f.campelo@bristol.ac.uk) - Principal investigator  
[Lindeberg Leite](mailto:lindpessoa@gmail.com) - First author
