# PoissonTimeVaryingInfo
Matlab code for calculating time varying information of single units spike patterns about stimulus ID and stimulus category
This code was developped for analysing electrophysioloical recordings of single auditory cortical neurons responding to playback
of vocalizations. This code can be applicable to any system where one wants to measure how much information a neuron convey in
its time-varying spike pattern about the identity of input stimuli or about classes or categories of stimuli.

The core functions of this repository are:

-> *info_poisson_model_calculus.m* which calculates the mutual information between n events (stimuli, categories of stimuli, motor output...)
and a poisson response defined by its n means.

-> *cumulative_info_poisson_model_calculus_MCJK.m* which calculates the cumulative information in a neural response given
conditional probability distributions of responses given stims at each time point. Although it is assumed that this stimulus
conditional probability is time independent (i.e know that stim identity and the rate at a previous time does not tell you 
about the rate at the current time - as in a poisson distribution), the joint probability distribution is not, making the
cumulative information calculation difficult. This algorithm uses a Monte Carlo approximation to estimate this joint probability
distribution. The error on the calculation of cumulative information is estimated by a Jack-knife procedure.

-> *info_cumulative_model_Calculus.m* which calulates the cumulative information in a neural response given
conditional probability distributions of responses given stims at each time point. Although it is assumed that this stimulus
conditional probability is time independent (i.e know that stim identity and the rate at a previous time does not tell you 
about the rate at the current time - as in a poisson distribution), the joint probability distribution is not, making the
cumulative information calculation difficult. This algorithm offers different methods to estimate this joint probability
distribution.
3 exact calculation methods that each have their own inconveniences: 'Exact_Mem' uses the computer memory and will crash
if the dataset is too large; 'Exact_Paths' do all the calculations one by one and will take infinite time if the dataset
is too large; 'Exact_HardDrive' will also perfom all the calculations, but saving temporary calculations on the hard drive,
eventualy crashing for large dataset.
2 approximation methods: 'MarkovChain' uses the parameters enterd as MarkovParameters to estimate the entropy of the neural
response with a history and a pace for the Markov chain given as inputs; 'MonteCarlo' uses a MonteCarlo estimation with a
weight correction, the number of samples used for the estimation is given an input.
