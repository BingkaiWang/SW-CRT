# R code for the paper "How to achieve model-robust inference in stepped wedge trials with model-based methods?"

I. The preprint paper is available at https://arxiv.org/abs/2401.15680. 

II. The folder "Simulation code" contains code for our simulation studies.

conSimulation: code for our simulation studies with continuous outcomes.

- List of Supporting Files: These supporting files are sourced in the main files that reproduce the numbers in the submitted manuscript and supporting information, for continuous outcomes.

1. SV.R = function to calculate the sandwich variance using proposed methods with a simple exchangeable correlation structure;
2. SV_NE.R = function to calculate the sandwich variance using proposed methods with a nested exchangeable correlation structure;
3. contGEE_BCV2.R = function to calculate the GEE-style bias-corrected sandwich variances with a simple exchangeable correlation structure (GEE-style bias-corrected sandwich variances were not evaluated in our paper, due to the focus on relatively large number of clusters);
4. contGEE_BCV2_ED.R = function to calculate the GEE-style uncorrected and bias-corrected sandwich variances with an exponential decay correlation structure (see Web Appendix D.3 for details; GEE-style bias-corrected sandwich variances were not evaluated in our paper, due to the focus on relatively large number of clusters).

- List of Main Files: These main files are used to reproduce the results in the submitted manuscript and supporting information, for continuous outcomes.

5. simuSWD_A1.R = reproduce results for Simulation Scenario A1 with a simple exchangeable or nested exchangeable correlation structure;
6. simuSWD_A1ed.R = reproduce results for Simulation Scenario A1 with an exponential decay correlation structure;
7. simuSWD_A2.R = reproduce results for Simulation Scenario A2 with a simple exchangeable or nested exchangeable correlation structure;
8. simuSWD_A2ed.R = reproduce results for Simulation Scenario A2 with an exponential decay correlation structure;
9. simuSWD_B1.R = reproduce results for Simulation Scenario B1 with a simple exchangeable or nested exchangeable correlation structure;
10. simuSWD_B1ed.R = reproduce results for Simulation Scenario B1 with an exponential decay correlation structure;
11. simuSWD_B2.R = reproduce results for Simulation Scenario B2 with a simple exchangeable or nested exchangeable correlation structure;
12. simuSWD_B2ed.R = reproduce results for Simulation Scenario B2 with an exponential decay correlation structure.

binSimulation: code for our simulation studies with binary outcomes.

- List of Supporting Files: These supporting files are sourced in the main files that reproduce the numbers in the submitted manuscript and supporting information, for binary outcomes.

1. sim_estimand_1.R = function to calculate estimands for Simulation Scenario C1;
2. sim_saturated_1.R = function to analyze data for Simulation Scenario C1 with a simple exchangeable or nested exchangeable correlation structure, without covariate adjustment;
3. sim_saturated_1x.R = function to analyze data for Simulation Scenario C1 with a simple exchangeable or nested exchangeable correlation structure, with covariate adjustment;
4. contGEE_BCV2_ED.R = function to calculate the GEE-style uncorrected and bias-corrected sandwich variances with an exponential decay correlation structure (see Web Appendix D.3 for details; GEE-style bias-corrected sandwich variances were not evaluated in our paper, due to the focus on relatively large number of clusters);
5. sim_saturated_1_ed.R = function to analyze data for Simulation Scenario C1 with an exponential decay correlation structure;
6. sim_estimand_2.R = function to calculate estimands for Simulation Scenario C2;
7. sim_saturated_2.R = function to analyze data for Simulation Scenario C2 with a simple exchangeable or nested exchangeable correlation structure, without covariate adjustment;
8. sim_saturated_2x.R = function to analyze data for Simulation Scenario C2 with a simple exchangeable or nested exchangeable correlation structure, with covariate adjustment;
9. sim_saturated_2_ed.R = function to analyze data for Simulation Scenario C2 with an exponential decay correlation structure.

- List of Main Files: These main files are used to reproduce the results in the submitted manuscript and supporting information, for binary outcomes.

10. simuSWD_C1.R = reproduce results for Simulation Scenario C1 with a simple exchangeable or nested exchangeable correlation structure, without covariate adjustment;
11. simuSWD_C1x.R = reproduce results for Simulation Scenario C1 with a simple exchangeable or nested exchangeable correlation structure, with covariate adjustment;
12. simuSWD_C1_ed.R = reproduce results for Simulation Scenario C1 with an exponential decay correlation structure;
13. simuSWD_C2.R = reproduce results for Simulation Scenario C2 with a simple exchangeable or nested exchangeable correlation structure, without covariate adjustment;
14. simuSWD_C2x.R = reproduce results for Simulation Scenario C2 with a simple exchangeable or nested exchangeable correlation structure, with covariate adjustment;
15. simuSWD_C2_ed.R = reproduce results for Simulation Scenario C2 with an exponential decay correlation structure.

III. The folder "data-analysis" contains code for reproducing our data application. The data we used are publicly available at https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/NSKFK2.
