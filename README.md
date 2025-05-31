# EpitopeTransfer

This repository contains data and code for the paper **EpitopeTransfer:
Improving linear B-cell epitope prediction through phylogeny-aware
transfer learning** by Lindeberg Pessoa Leite, Teófilo E. de Campos,
Francisco Pereira Lobo and Felipe Campelo.

## Table of Contents

-   [Dependencies](#Dependencies)
-   [Running models for published metrics](#Running-models-for-published-metrics)
-   [Running the analyses](#Running-the-analyses)
-   [Contact](#contact)

## Dependencies 

This project requires Python 3.10, although other versions are likely
compatible. Follow the instructions below to set up your environment.

1.  **Install Python 3.10**

Ensure that Python 3.10 is installed on your system. If not, you can install it using your system's package manager. For Ubuntu-based distributions, you can use:

``` bash
    sudo apt update
    sudo apt install software-properties-common -y
    sudo add-apt-repository ppa:deadsnakes/ppa -y
    sudo apt update
    sudo apt install python3.10 python3.10-venv python3.10-distutils -y
```

2.  **Activate epitopetransfer environment**

``` bash
  source esm2/bin/activate
  # or
  source esm1b_v1/bin/activate
  # or
  source esm1b_v2/bin/activate
```
To exit any environment, run: deactivate

## Running models for published metrics

To run the models for specific taxa or for all taxa included in the
study, use the following command format in the terminal:

``` bash
  python3.10 main.py [base_model] [taxa]
```

Replace [base_model] to esm1b or esm2, and [taxa] to the desired taxa or all (for esm2 only). Example:

``` bash
  source esm2/bin/activate
  python3.10 main.py esm2 all # (the 'all' option is available for esm2 base model only)
```
To exit the `esm2` environment, use this command: `deactivate`

``` bash
  source esm1b_v1/bin/activate
  python3.10 main.py esm1b bpertussis
```

``` bash
  source esm1b_v2/bin/activate
  python3.10 main.py esm1b mononegavirales
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


**Note:**  For the esm1b base model, activate the esm1b_v1 environment (source esm1b_v1/bin/activate) for taxa in the first column, or the esm1b_v2 environment (source esm1b_v2/bin/activate) for taxa in the second column. For the esm2 base model, activate the esm2 environment (source esm2/bin/activate) for all taxa. 

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
