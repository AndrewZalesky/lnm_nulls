# Permutation-based null models for lesion network mapping (LNM)

[Lesion network mapping (LNM)](https://academic.oup.com/brain/article-abstract/138/10/3061/2468715) is a popular methodology to to identify brain circuits disrupted by spatially distributed lesions associated with the same symptom or clinical phenotype.

We evaluate the statistical power of spatial and topological null models for the simplified LNM model proposed by [van den Heuvel and colleagues (2026)](https://www.nature.com/articles/s41593-025-02196-7).

## Simplified LNM model
The simplified LNM is based on a single functional connectivity matrix defined at a regional scale. There is no modeling of individuals. Each lesion maps to one and only one region. Depending on the lesion assignment type, multiple lesions can be mapped to the same region. This simplified LNM formulation is best viewed as an illustrative model and may differ from the properties of voxel‐resolution LNM applied across individuals.

## Simulated functional connectivity matrix
The matrix is generated using a degree-controlled stochastic block model (dcSBM). Connectivity weights are sampled from a beta distribution. The matrix is symmetric, undirected, fully connected and comprises *B* modules. Each module contains *M* nodes (regions) and the total number of nodes is thus *N=BM*. The node strength distribution is approximately log-normal.  

## Lesions sets and ground truth lesion network map
The ground truth lesion network map is always defined by a single module comprising the simulated functional connectivity matrix. Rather than conforming to a single predefined module at one scale, lesion network maps in practice may comprise more complex or multiscale topological organization. We consider a total of *K* lesions. A proportion *alpha* of these *K* lesions is necessarily assigned to nodes comprising the ground truth module. The remaining lesions are uniformly distributed at random across all nodes, irrespective of modular allegiance. The parameter *alpha* determines the propotion of lesions in the lesion network map and we plot statistical power as a function of *alpha*. 

- When *alpha=0* all lesions are distributed uniformly and this case evaluate the specifcity (i.e. family-wise error rate, FWER) of the methodology.
- As *alpha* increases, more lesions are concentrated in the ground truth LNM and statistical power thus increases.

## Null models 
We evaluate wpatial and topological null modelw. The spatial null model is based on the null hypothesis that the observed lesion network map is consistent with an equal number of equally sized lesions distributed uniformly at random throughout the brain. The topological null hypothesis is that the observed lesion network map is consistent with a functional network in which higher-order topological properties hypothesized to support disease-specific circuits are absent. To generate samples from the topological null model, we draw networks from the same dcSBM used to generate the functional connectivity matrix, while explicitly excluding modular structure from the generative process. The null networks are matched in strength distribution to the simulated functional connectivity matrix. 

## Statistical testing 
The null distribution for both null models is constructed by identifying the node with the maximum value in the lesion network map and storing this value. This is a conservative null distribution that ensures strong control of the FWER under assumptions of exchangeability. Nodes in the observed lesion network map are deemed statistically significant if they exceed the 95th percentile of the null distribution, ensuring control of the FWER at 0.05. The number of permutations and samples in the null distribution is determined by *Perms*.  

## Evaluation metrics
We generate many lesions sets and corresponding lesion network maps for a given set of parameters (i.e. *N*, *alpha*, *K*, etc) and undertake the null hypothesis testing described above. Each such repeated iteration is referred to as a trial. The following measures are computed as the proportion of certain events across trials.  

- **Statistical power** is the proportion of trials for which one or more nodes are deemed significant with the ground truth lesion network map when *alpha>0*.
- **Family-wise error rate (FWER)** is the proportion of trials for which one or more nodes are deemed significant anywhere in network when *alpha=0*.

## Matlab code

The Matlab code provided here enables implementation of the simulations to generate components of the below figure.
- **main.m** is a script to evaluate all metrics as a function of *alpha*, *K* or a single parameter combination. The choice of null model, network parameters as well as the number of permutations and trials is specifcied in this script. 
- **dcsbm.m** generates weighted, symmetric undirected networks with and without modular structure.
- **lesion_assignment.m** maps lesions to nodes based on module allegiance. 
- **lnm_compute.m** computes the lesion network map and contains the for-loop to generate samples to build the null distribution.
- **precompute_topology_randomization.m** generates a set of functional connectivity matrices under the null hypothesis. The same set of precomputed matrices are used for all trials to save computation time.
- **show_progress.m** is a helper function.

Note that the [cbrewer package](https://github.com/scottclowe/cbrewer2) is needed to generate some figures. 

<br><br>


![figure](https://github.com/user-attachments/assets/0fbcbe6d-2b32-45de-92df-787b696296e9)
