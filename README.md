# CircadianBarley

This folder contains the codes to run system identification between genes from gene expression time series data
and perform the nu-gap analysis which assesses the consistency of the models in different experimental conditions (mutation).

The "main.m" function provides examples on 'How to' run the network inference. 
For the details, please refer to the Supplemental Information of the original submission.

1) [fitness, models, IO, AIK, dcs] = just_tfest(orders, Ts, data);

This program estimates models using every possible pairwise combination in the 'data' provided as input-output sets. The order of the models is specified by the variable 'orders'. Several orders can be chosen at the same time as following: [1 2].
The variable 'data' is supposed to be given in transcripts/time points (row/colomn). The sampling frequency 'Ts' is 1 by default and does not influence the fitness score.

The code calls 'tfest' and 'compare' from the System Identification Toolbox to respectively compute the models and the fitness scores.

The outputs are 
'fitness': a matrix with the fitness of model (i,j) using data in row j as input and data in row i as output
'models': a matrix with the transfer functions (i,j) that are derived from the models as above
'AIK: a matrix with the information criterion (Aikake) of the models.
'dcs': a matrix with the DC gain of the models. The latter can be used to assess the sign of the regulation.

2) [nugap,freq,w] = richgap(L1, L2, p1, p2)

This program computes the nu-gap between two transfer functions (L1 and L2), considering only frequencies between p1 and p2.
For instance, if the interesting periods lie between 18 and 36, then the correct p1 and p2 inputs should be p1=2*pi/18 and p2=2*pi/36

The outputs of the code are :
'nugap': a close accurate of the value of the nu-gap between the transfer functions
'freq': the frequency to which the nugap correpond (the frequency that correspond to the maximum distance between the two transfer functions in the Riemann plane)
'w': the set of frequencies considered when computing the nu-gap estimation

NOTE: the code does not check if the winding number condition is satisfied. Use the MATLAB(TM) code '[gap,nugap] = gapmetric(p0,p1)' for that. 
