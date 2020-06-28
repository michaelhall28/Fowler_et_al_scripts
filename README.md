Python scripts and notebooks for Fowler et al.  

## Requirements  
All scripts and notebooks require Python 3.  

The notebooks (.ipynb) analyse mutations in a subset of proteins.   
These require the installation jupyter python package (https://jupyter.readthedocs.io/) and the code from github repository https://github.com/michaelhall28/darwinian_shift.

To run the comparison with Clinvar annotated variants, download and uncompress the variants_summary.txt file from Clinvar (https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz - accessed April 2020).

The python script competition\_simulation\_script.py runs stochastic simulations to see how much inter-sample variation in observed mutation counts might be expected under uniform environmental conditions.
This script requires the installation of the code from github repository https://github.com/michaelhall28/clone-competition-simulation   

