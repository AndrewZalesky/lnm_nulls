# Permutation-based null models for clinical lesion maps 

Here we apply our spatial null model to the clinical lesions sets [made available by van den Heuvel et al](https://github.com/dutchconnectomelab/lesionnetworkmapping/tree/main)

Our spatial null model approximately preserve regional strength (i.e., node degree) when randomizing lesion loci. 

- *main.m* is script that implements both the compact LNM formulation as well as permutation testing

Note that the functional connectivity matrix is first corrected for know relationship between functional connectivity and distance. 
