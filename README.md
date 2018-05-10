# Ebola transmission dynamics in the DRC

# File overview
## Code
 - This folder contains all of the code necessary to run the analyses for the paper... For the analysis we fit two different epidemiological models (intervention and superspreading) to multiple DRC outbreaks. We estimate the parameters from the models, and then
## Data

# Running the code

This script uses the following scripts:
1. 
2. 

Questions:
1. Can we abstract process from the mif2_run function?
2. Why are there two mif runs in that function?


Comments
1. Overall need all like-files to be easy to understand their connection
2. Also need all file connections to be clear and more transparent
  - e.g. main_functions script should be able to be used by each mod we throw at it, if it is specific to the model we use, should be in the specific model script and associated
3. Want functions to be clear and transparent
4. One function - one task
5. What is the second liks_global step of mif2 for?

6. Probably stop it from using the clusters?




2 Different models
6 Different outbreaks

12 different combinations

Figure 1 - comparison between two models for each outbreak - showing estimates and fits

- To get this figure what do we need?
- mif fits for each model/outbreak combination
- confidence intervals for each model/outbreak/parameter combination


Fig 2 - TBD