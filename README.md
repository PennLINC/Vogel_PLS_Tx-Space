# Vogel_PLS_Tx-Space
Data and code for Vogel, Seidlitz et al., Transcriptomic Gradients paper

Included are several Jupyter notebooks running through each of the analyses in the paper, showing the exact code and data used to prepare and run each analysis and sub-analysis, and generate the figures from the manuscript.

The Notebooks are divided thematically, such that each corresponds to a certain set of analyses. Due to the size of some of the files, you will need to run some of these notebooks (particularly NB1, NB2 and NB3) in order to generate the data used in later notebooks. I tried to demarcate areas where a previous notebook must be run in order to complete an analysis. In addition, some data will needed to be downloaded from the public repositories on the internet. Instructions and links are included in the notebook. For best results, you should place all downloaded data into the `data/` folder of this repository.

# Requirements
Running these notebooks will require Python 3.7.3 or above, Jupyter, and a number of Python packages. A `requirements.txt` file containing all of the necessary elements is included in the git repo. 

# Reproducing analyses on your machine
To reproduce these analyses on your computer, all you have to do is:

* Install python if you haven't aready.
* Clone or download this entire repository.
* Use the requirements.txt file to recreate my python environment with the following commands (inside this git repo):

```
virtualenv pls_gxp
source plx_gxp/bin/activate
pip install -r requirements.txt
```

You must be sure to be "inside" the pls_gxp environment to run these notebooks. You can enter the environment at any time by running `source plx_gxp/bin/activate` from within the git repo, and you can leave the environment at any time by typing `deactivate`. You'll know you're in the environment if you see `(pls_gxp)` in the very front of your command prompt.

To run a notebook, navigate to repo and run the notebooks: jupyter notebook <Notebook.ipynb>

# Troubleshooting
If there are any problems, please don't hesitate to raise an issue.
