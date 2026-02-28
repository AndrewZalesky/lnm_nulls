# Permutation-based null models for clinical lesion maps 

Here we apply our spatial null model to the clinical lesions sets [made available by van den Heuvel et al](https://github.com/dutchconnectomelab/lesionnetworkmapping/tree/main)

Our spatial null model approximately preserves regional strength (i.e., node degree) when randomizing lesion loci. 

- *main.m* is a script that implements both the compact LNM formulation as well as permutation testing
- *compute_lnm.m* computes the lesion network map
- *compute_null.m* compute null computes one sample from the null distribution
- All other functions are helpers. 

Note that the functional connectivity matrix is first corrected for the known relationship between functional connectivity and distance using regression. 
