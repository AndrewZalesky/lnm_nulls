# Permutation-based Null models for lesion network mapping (LNM)

Lesion network mapping (LNM) is a popular methodology to to identify brain circuits disrupted by spatially distributed lesions associated with the same symptom or clinical phenotype.

We evaluate the statistical power of spatial and topological null models for the simplified LNM model proposed by van den Heuvel and colleagues (2026).

##Simplified LNM model
The simplified LNM is based on a single functional connectivity matrix defined at a regional scale. There is no modeling of individuals. Each lesion maps to one and only one region. Depending on the lesion assignment type, multiple lesions can be mapped to the same region. This simplified LNM formulation is best viewed as an illustrative model and may differ from the properties of voxel‐resolution LNM applied across individuals.

##Simulated functional connectivity matrix
The matrix is generated using a degree-controlled stochastic block model (SBM). Connectivity weights are sampled from a beta distribution. The matrix is symmetric, undirected, fully connected and comprises *B* modules. Each module contains *M* nodes (regions) and the total number of nodes is thus *N=BM*. The node strength distribution is approximately log-normal.  

##Lesions sets and ground truth lesion network map
The ground truth lesion network map is always defined by a single module comprising the simulated functional connectivity matrix. Rather than conforming to a single predefined module at one scale, lesion network maps in practice may comprise more complex or multiscale topological organization. We consider a total of *K* lesions. A proportion *alpha* of these *K* lesions is necessarily assigned to nodes comprising the ground truth module. The remaining lesions are uniformly distributed at random across all nodes, irrespective of modular allegiance. The parameter *alpha* determines the propotion of lesions in the lesion network map and we plot statistical power as a function of *alpha*. 

- When *alpha=0* all lesions are distributed uniformly and this case evaluate the specifcity (i.e. family-wise error rate, FWER) of the methodology.
- As *alpha* increases, more lesions are concentrated in the ground truth LNM and statistical power thus increases.

##Null models 
We evaluate a spatial and topological null model. 




The Matlab code provided here enables implementation of the simulations to generate components of the below figure. 



![figure](https://github.com/user-attachments/assets/0fbcbe6d-2b32-45de-92df-787b696296e9)
