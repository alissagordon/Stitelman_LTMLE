# Stitelman_LTMLE
Toy model implementing the LTMLE methodology in Stitelman, et al (2012) Doi: 10.1515/1557-4679.1334

The current LTMLE package employs the method of iterative conditional expectations. However, an alternative approach for TMLE of longitudinal data exists in this paper by instead targeting the entire likelihood of the data. This approach has the estimation benefits over the ICE approach, but code that implements it is not readily available. The purpose of this repository is to demonstrate how one may employ the likelihood approach of LTMLE. While we only code this approach for three time points (t=3), it is a good starting point for those who may be interested in using this methodology. 

Our approach and simulation results are available in the slides. While we do not employ the Markov assumption, it is heavily recommended under more timepoints and complex time-varying covariates. See the paper for more detailed explanation of the approach.

Three files are attached. First, there are the slides that explain our approach and show simulation results after comparing this approach against other estimation techniques. Second, there is R code showing a single implementation of the Stitelman LTMLE approach. Finally, there is an ipynb file that was used for the simulations in the slides. The ipynb file was coded to use R.

Authors of the code and slides in the repo: Alissa Gordon, Yilong Hou, Kaitlyn Lee, and Sylvia Song 
Affiliation with UC Berkeley Biostatistics Department

For questions and inquiries: alissa_gordon@berkeley.edu
