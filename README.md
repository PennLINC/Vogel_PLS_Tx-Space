# Vogel_PLS_Tx-Space
Data and code for [Vogel et al. 2024, Transcriptomic Gradients paper](https://www.pnas.org/doi/10.1073/pnas.2219137121), now accepted at *PNAS*. See [project page](https://pennlinc.github.io/Vogel_PLS_Tx-Space/) for additional information.

Included are several Jupyter notebooks running through each of the analyses in the paper, showing the exact code and data used to prepare and run each analysis and sub-analysis, and generate the figures from the manuscript.

The Notebooks are divided thematically, such that each corresponds to a certain set of analyses. Due to the size of some of the files, you will need to run some of these notebooks (particularly NB1, NB2 and NB3) in order to generate the data used in later notebooks. I tried to demarcate areas where a previous notebook must be run in order to complete an analysis. In addition, some data will needed to be downloaded from the public repositories on the internet. Instructions and links are included in the notebook. For best results, you should place all downloaded data into the `data/` folder of this repository.

The following lists which notebooks are needed to reproduce analyses displayed in each manuscript figure:

NB01 --> None (data generation)

NB02 --> Fig 1B,E,F,G,H; Fig S1A,C,D,F; Fig S3

NB03 --> Fig 1C,D; Fig S1E; Fig S4 (except panel E), Table S1, Table S7

NB04 --> Most of Fig 2; Fig S4E; Fig S6; Fig S6; Fig S8

NB04a --> The rest of Fig 2

NB05 --> Fig 1K, Table S2

NB06 --> Fig 1J; Fig 4; Fig S11

NB07 --> Fig 4C; Fig S11A; Fig S12

NB08 --> Fig 1L; Fig 3; Fig S9

NB09 --> Fig 5, Fig S13, Table S3, Table S4, Table S5, Table S6

# Requirements
Running these notebooks will require Python 3.7.3 or above, Jupyter, and a number of Python packages. A `requirements.txt` file containing all of the necessary elements is included in the git repo. 

# Reproducing analyses on your machine
To reproduce these analyses on your computer, all you have to do is:

* Install python if you haven't aready. We suggest Anaconda.
* Clone or download this entire repository.
* Use the requirements.txt file to recreate my python environment with the following commands (inside this git repo):

```
conda create -n pls_gxp python=3.7.3
conda activate pls_gxp
pip install -r requirements.txt
```

You must be sure to be "inside" the pls_gxp environment to run these notebooks. You can enter the environment at any time by running `source plx_gxp/bin/activate` from within the git repo, and you can leave the environment at any time by typing `deactivate`. You'll know you're in the environment if you see `(pls_gxp)` in the very front of your command prompt.

To run a notebook, navigate to repo and run the notebooks: jupyter notebook <Notebook.ipynb>

# Troubleshooting
If there are any problems, please don't hesitate to raise an issue.
