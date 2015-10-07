# IKAP (Inference of kinase activities from Phosphoproteomics)
Matlab algorithm to estimate kinase activities from Phosphoproteomics data

Step 1-1 SearchMotifs:
The first function searches the data for motifs that are known from the PSP database to be phosphorylated by a certain kinase. The database can be downloaded from http://www.phosphosite.org/staticDownloads.do. It should be reduced to one column containing the kinase names and one containing the motifs. Since the exact position within the motif is of great importance, the function compares the strings which are left and right from the central amino acid to those of the motifs in the database. It adds a column with the kinases found to the data matrix and takes as an input the original dataset as a cell array, the PSP database as a two-column cell array and the column containing the motif in the dataset. The data matrix should begin with a column containing the protein IDs (e.g. gene names) followed by a column for the motifs and a variable number of data columns. Thhe phosphorylated amino acid should be in the middle of the motif.

Step 1-2 MakeKin:
This function creates a list of the kinases that are known to phosphorylate at least one motif in the dataset. If additionally proteome data are available (optional 3rd argument, list of protein names) the function reduces the kinase list to those which are present in the proteome. Finally, it produces a new data matrix which only contains those phosphosites for which a kinase is known.

Step 1-3 CreateTT:
The third function creates a truth table out of the reduced dataset, with 1 for true interactions and 0 for false. This will be used as a mapping in the next step to estimate the kinase activities and affinity parameters on the basis of only those measurement values, the respective kinases have an influence on. 

Step 1-4 FitActivities (using ComputeCostGlobal, ComputeCostLocal and APReadout):
FitActivities estimates the activities of the kinases in kin using the Matlab built-in function fmincon. The number of iterations can be specified using the input argument iter. Each iteration starts with randomized parameter vectors and cycles between an optimization of the affinity parameters by means of the global cost function and an optimization of the kinase activities for each condition by means of the local cost function. Since the gradients for the parameters are delivered by the cost functions, the ‘GradObj’ argument in the option settings should be turned to ‘on’. The while-loop continues as long as the difference between both costs is larger than a certain value, in this case 10. The upper and lower bounds for parameter estimation as well as the starting values should be modified according to the user´s needs. The output variables AP and K contain the estimated affinity parameters and kinase activities, cost is a list of all calculated costs with mincost being the lowest. This function can also be executed on multiple cores using the parallel computing toolbox. In this case, the first for-loop needs to be replaced by a parfor-loop.

The global cost function calculates the cost over all conditions and the gradients for the affinity parameters. The global cost is defined as the sum of the squared distances of the measured phosphosites and the mean of the products of activity and affinity of all kinases phosphorylating the phosphosite, over all conditions. The gradient for each affinity parameter is given as its partial derivative of the cost function.

The local cost function calculates the cost for one condition and derives the gradients for the kinase activities within the respective condition. As in the global case, the cost is defined as the sum of the squared distances of the measured values and the means of the activity*affinity products of the assigned kinases, but without summing over all conditions.

This function creates a cell array that allocates the estimated affinity parameters to the respective kinase-target-links.


Step 1-5 IdentKin (using FitIdent) (optional):
IdentKin.m performs an identifiability analysis of the estimated kinase activities using the fitting function FitIdent.m, which makes use of fmincon. The activities at the condition c of the kinases given in ks (numbers as they are listed in kin) are decreased and increased in a step-wise manner (increment ranging from -3 to 3) and fixed, followed by an optimization of the remaining parameters (both affinities and activities). The function outputs the matrices plJ, containing the obtained costs, and plk, containing the activity values with which the optimizations were performed. By applying a parfor-loop, it can be executed in parallel.

This function performs the optimization task within the identifiability analysis.

Step 2-1 ComputeDistances:
The first function of part 2 computes correlation coefficient based p-values for all kinase-phosphosite combinations. The resulting distance matrix is a cell array with kinases on the x-axis and phosphosites on the y-axis. In addition it produces the matrix valids, which contains the number of valid measurement values that went into the calculation of each p-value. 

Step 2-1-2 MakeDataVal (optional):
This function generates a dataset from which the links that have an annotated kinase in PSP are removed. This should be done previous to validation.

Step 2-2 MakePsig:
Psig is a list of all kinase-substrate pairs that have a p-value below the given false discovery rate q, calculated with the Benjamini-Hochberg method.

Step 3-1 MakePnsig:
This function is used to create a list of the same size as psig but with randomly selected links. By comparing psig and pnsig within the following functions, it is possible to validate the newly found kinase-substrate links.

Step 3-2 ValidDB:
This function searches the database of your choice for the kinase-target-links in psig and pnsig and calculates a p-value for the positive findings based on a fisher exact test. If the pair is found, the function assigns a 1 to the respective row in the fifth column of psig or pnsig. It can be executed in parallel. If no parallelization is used, the parfor-loop needs to be replaced by a normal for-loop.

Step 3-3 ValidMOT
This function searches the NetworKIN database for the motifs in psig and pnsig. If the motif is found, it notes the respective likelihood in the sixth column of psig or pnsig. Based on these columns, the function calculates the p-values pSUM, which is based on the sum of the likelihoods in psig and pnsig, and pNUM, which corresponds to the number of likelihoods that are larger than 1. The first loop can be executed in parallel.

