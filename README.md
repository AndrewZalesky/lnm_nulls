# Permutation-based Null models for lesion network mapping (LNM)

Lesion network mapping (LNM) is a popular methodology to to identify brain circuits disrupted by spatially distributed lesions associated with the same symptom or clinical phenotype.

We evaluate the statistical power of spatial and topological null models for the simplified LNM model proposed by van den Heuvel and colleagues (2026). The simplified LNM is based on a single functional connectivity matrix defined at a regional scale. The matrix is generated using a degree-controlled stochastic block model (SBM). Connectivity weights are sampled from a beta distribution. The matrix is symmetric, undirected, fully connected and comprises *B* modules. Each module contains *M* nodes (regions) and the total number of nodes is this *N=BM*.  

The Matlab code provided here enables implementation of the simulations to generate components of the below figure. 



![figure](https://github.com/user-attachments/assets/0fbcbe6d-2b32-45de-92df-787b696296e9)
