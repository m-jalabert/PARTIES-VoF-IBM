# Documentation for Oscillation-Testcases

The oscillation testcases are divided into the oscillation of the container/domain and the oscillation of the particle. 


## Oscillation of the domain

This testcase represents a fluid container oscillating unidirectional in zero gravity according a cosine-function.
A single particle is initially located to the center of the domain and oscillates due to the fluid-particle interaction. 
The setup as well as the fluid and particle characteristics are chosen based on the experiments of LÉspérance et al. (2005). 
Those experiments were chosen to validate the ratio of the particle to fluid excursion, as shown in Fig. 3 of the article.
The following spreadsheet shows the dimensional values of the experiments, the non-dimensional values of the validation, and the simplified non-dimensional values of the present testcase:

 Variables | Experiments | Validation Setup [-] | Testcase Setup [-] 
:-: | :-: | :-: | :-: |
Dp | 3.96E-3 [m] | 1 | 1 |
Lx / Ly / Lz | 3.96E-2 [m]  | 10 | 10 / 10 / 3
cells / Dp | - | 20 | 10
rho_p / rho_f | 4.68 | 4.68 | 4.68 
nu_f | 4.00E-5 [m2/s] | 3.64E-2 | 3.64E-2
frequency | 70 [Hz] | 1 | 1
amplitude | 3.96E-4 [m] | 0.1 | 0.1

For the dimensionless values, the particle diameter Dp [m] was chosen as length scale, the frequency f [Hz] as time scale, and the fluid density rho_f [kg/m3] as density scale.
Based on the non-dimensional values of the validation, the dimensionless size of Lz was reduced from 10 to 3 and the grid resolution was reduced from 20 cells/Dp to 10 for the present testcase.


## Oscillation of the particle

In this testcase, the fluid container is fixed and the particle oscillates in x-direction according a prescribed sine-function. 
The same simplified setup is chosen as shown in the spreadsheet above. 

