# cell_lineage_reconstruction_using_felsensteins_and_MCMC
Cell-Lineage Reconstruction using Felsenstein's and MCMC

FelsenMCMC_Main is the script that runs the algorithm calling the functions felsensteintree and allposparents

RF_scoring is an r script that will automatically score the output results using the function RF.dist from the phangorn package

allBarcodes is used to initialize the transition matrix in FelsenMCMC_Main

train_setDREAM2019 contains the ground truth trees in newick format

This work was completed with Zhonghao Dai and Weida Ma and the final report is attached for those interested in our project. 

The idea stemmed from a challenge proposed by the Allen Institute in partnership with DREAM, and the data used was originally from the Elowitz Lab at Caltech. (Link: https://www.synapse.org/#!Synapse:syn20692755/wiki/595096) 
