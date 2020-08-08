# PLGA-Model

Collaborated with Dr. Anil Shrirao

###### Description ######

This is a devised mathematical MATLAB model simulating the drug release profile of lidocaine from PLGA nanoparticles over time.
Based on the mechanistic properties of drug release, the properties of drug diffusion and degradation are taken into account
with the simulation of drug release from the nanoparticle construct. This can be optimized in the future to fit different 
polymers and loaded drugs by changing chemical and physical properties of the desired polymer/drug of choice.


##### Folder Contents ###### 

PLGAModel.mlapp:

This is the GUI application for the model. This maintains all the front-end control of everything that is going on. By
running the app and entering relevant information, the information output from the model is displayed both as a plot 
and as a data table. The simulation takes approximately ~5 minutes.


Degrado_main1.m:

This file is the main script containing key information such as constants, coefficients, and other physical/chemical 
properties pertaining to PLGA and the properties of the drug. Furthermore, it takes into account initial conditions
for all 14 differential equations and data manipulation of solution obtained from solving the system of equations.
Mass fluxes and average values are also calculated in this script.


Degrado_main2.m:

This file is exactly the same as Degrado_main1 but is set up as a function as opposed to a script. This function
takes in loading capacity and drug molecular weight as parameters. This file is primarily used for GUI compatibility.


degrado_ode1.m:

The purpose of this file is to solve the set of differential equations. This takes into account all the equations of
diffusion and degradation described in the paper including their associated boundary conditions. Intermediate derivative
calculations are done as well that are relevant to the simplified expression of the set of degradation equations.


fitting_Ddrug1.m:

This file serves the purpose of calculating the diffusion coefficient of the drug using the process of solving 
nonlinear least-squares problems. This involves an iterative process of data fitting through solving a system 
of ordinary differential equations multiple times (identical setup as that in degrado_ode1.m) where the stopping
condition is when the difference between the values of the experimental values and the values of the model are 
below a set tolerance level. This outputs the diffusion coefficient of the drug.


Siepmann_et_al_exp_data.xls:

This Excel file holds the experimental data from the diffusion study. This file is read and loaded into MATLAB.


Siepmann_et_al_sim_drug_release.xls:

This Excel file holds the output data from running the MATLAB script. Data is automatically written to the file
once the script is run.
