function [tempo,Rel_sim,time_exp,Rel_exp] = Degrado_main2(PMdrug,Drug_load)

% Degrado main
clc 
clear all

global Radius Kd Radius_Fragments deltaR Sh Dmon_pol Dmon_wat Dwat_wat Dwat_pol Cw_bulk options
global time IC PMmon pKa Cmin PMn beta time_exp Rel_exp
global Ddrug_wat Ddrug_pol

% Particle geometry
Diameter = 53*2*10^-4;              % Particle diameter, cm
Radius = Diameter/2;                % Paricle radius
SurfaceArea = 4*pi*Radius^2;                % Surface
volume = 4/3*pi*Radius^3;           % Volume

% Physical and chemical data
% Polymer Table1
PMw = 32440;                        % TableS5.2 30*1000 Weight average molecular weight, g mol-1
PD = 1.6;                           % Polydispersity index, [-]
PMn = PMw/PD;                       % Number average molecular weight, g mol-1
rho_pol = 1.2;                      % Polymer density, g cm-3
PM_lact = 90.08;                    % Lactic acid molecular weight, g mol-1
PM_glyc = 76.05;                    % Glycolic acid molecular weight, g mol-1
PMmon = (PM_lact + PM_glyc)/2;      % Pseudomonomer molecular weight, g mol-1
Massa = rho_pol*volume;             % Particle mass, g

% Diffusion coefficients
Dmon_wat = 10^-5;                   % Diffusivity of monomer in water, cm2 s-1
Dwat_wat = 2.3*10^-5;               % Diffusivity of water in water, cm2 s-1
Dmon_pol = 10^-10;                  % Diffusivity of monomer in polymer, cm2 s-1
Dwat_pol = 10^-8;                   % Diffusivity of water in polymer, cm2 s-1
Ddrug_wat = 10^-6;                  % Diffusivity of drug in water, cm2 s-1
beta = 2.5;

% Kinetic constants
Kd = 0.00243170660852174;                     % Hydrolysis kinetic constant

% Drug Data
PMdrug = 234.34;                            % Drug molecular weight, g mol-1
Drug_load = 4;                              % Drug loading, % w/w
Drug_mass = rho_pol*Drug_load/100;          % Drug concentration, g cm-3 

% Sherwood number
Sh = 2;

% Water
Cw_bulk = 1/18;                     % Water molar density, mol cm-3

%  Acid dissociation Con stant
pKa = 3.83;

% Experimental data
time_exp = table2array(readtable('Siepmann_et_al_exp_data.xls', 'Sheet', 'Drug_release', 'Range', 'K3:K13')).*86400;
Rel_exp = table2array(readtable('Siepmann_et_al_exp_data.xls', 'Sheet', 'Drug_release', 'Range', 'L3:L13'))./100;

% Parameters for numerical integration
Radius_Fragments = 50;                        % Number of discretization points for radius discretization
deltaR = Radius/Radius_Fragments;             % Distance between grid points
Simulation_Days = 35;                        % Time to be simulated, days
Simulation_Time_Seconds = Simulation_Days*86400;             % Time to be simulated, seconds
time = 0:3600:Simulation_Time_Seconds;              % Output for time

% Initial conditions need to mentioned for all the interval for each
% depnedant variable
Cmin = 0;
Cm0 = zeros(Radius_Fragments,1);
Cw0 = zeros(Radius_Fragments,1);
mu00 = zeros(Radius_Fragments,1);
mu10 = zeros(Radius_Fragments,1);
mu20 = zeros(Radius_Fragments,1);
Cds0 = zeros(Radius_Fragments,1);

for i = 1:1:Radius_Fragments
   Cm0(i) = Cmin;                   %add your value by changing Cmin
   Cw0(i) = 0;                      %is menioned as 0 in paper
   mu00(i) = rho_pol/PMn;           % equation 11.a supplimentary information
   mu10(i) = PMn*mu00(i)/PMmon;     % equation 11.b supplimentary information
   mu20(i) = PD*mu10(i)^2/mu00(i);  % equation 11.c supplimentary information
   Cds0(i) = Drug_mass/PMdrug;      % 7.b standard equation of concentration=Drug Mass/Drug molecular weight
end

IC = [Cm0 Cm0 Cm0 Cm0 Cm0 Cm0 Cm0 Cm0 Cm0 Cw0 mu00 mu10 mu20 Cds0];
options = odeset('MaxStep',3600,'Stats','off','NonNegative',[1:1:14*Radius_Fragments]);

% Command for fitting solute diffusion coefficient
Ddrug_pol = fitting_Ddrug1();

% Numerical integration
[tempo ODE_Result_i] = ode15s(@degrado_ode1, time, IC, options);

% Retrieving the specific variables from the resulting matrix
Cm = zeros(length(tempo),Radius_Fragments);
C2 = zeros(length(tempo),Radius_Fragments);
C3 = zeros(length(tempo),Radius_Fragments);
C4 = zeros(length(tempo),Radius_Fragments);
C5 = zeros(length(tempo),Radius_Fragments);
C6 = zeros(length(tempo),Radius_Fragments);
C7 = zeros(length(tempo),Radius_Fragments);
C8 = zeros(length(tempo),Radius_Fragments);
C9 = zeros(length(tempo),Radius_Fragments);
Cw = zeros(length(tempo),Radius_Fragments);
mu0 = zeros(length(tempo),Radius_Fragments);
mu1 = zeros(length(tempo),Radius_Fragments);
mu2 = zeros(length(tempo),Radius_Fragments);
Cds = zeros(length(tempo),Radius_Fragments);

lambda = zeros(length(tempo),Radius_Fragments);
PMnum = zeros(length(tempo),Radius_Fragments);
PMw = zeros(length(tempo),Radius_Fragments);

for i = 1:1:length(tempo)
    for j = 1:1:Radius_Fragments
        Cm(i,j) = ODE_Result_i(i,j);
        C2(i,j) = ODE_Result_i(i,j + Radius_Fragments);
        C3(i,j) = ODE_Result_i(i,j + 2*Radius_Fragments);
        C4(i,j) = ODE_Result_i(i,j + 3*Radius_Fragments);
        C5(i,j) = ODE_Result_i(i,j + 4*Radius_Fragments);
        C6(i,j) = ODE_Result_i(i,j + 5*Radius_Fragments);
        C7(i,j) = ODE_Result_i(i,j + 6*Radius_Fragments);
        C8(i,j) = ODE_Result_i(i,j + 7*Radius_Fragments);
        C9(i,j) = ODE_Result_i(i,j + 8*Radius_Fragments);
        Cw(i,j) = ODE_Result_i(i,j + 9*Radius_Fragments);
        mu0(i,j) = ODE_Result_i(i,j + 10*Radius_Fragments);
        mu1(i,j) = ODE_Result_i(i,j + 11*Radius_Fragments);
        mu2(i,j) = ODE_Result_i(i,j + 12*Radius_Fragments);
        Cds(i,j) = ODE_Result_i(i,j + 13*Radius_Fragments);
    end
end

% Number and weight average molecular weights
PMnum = mu1./mu0.*PMmon; %equation 11.a-c
PMw = mu2./mu1.*PMmon;

% Calculation of average quantities as a function of time through integral average over radius
R = deltaR:deltaR:Radius;
Cmave = zeros(length(tempo),1);
C2ave = zeros(length(tempo),1);
C3ave = zeros(length(tempo),1);
C4ave = zeros(length(tempo),1);
C5ave = zeros(length(tempo),1);
C6ave = zeros(length(tempo),1);
C7ave = zeros(length(tempo),1);
C8ave = zeros(length(tempo),1);
C9ave = zeros(length(tempo),1);
Cwave = zeros(length(tempo),1);
Cdsave = zeros(length(tempo),1);
Cdstot = zeros(length(tempo),1);
mu0ave = zeros(length(tempo),1);
mu1ave = zeros(length(tempo),1);
mu2ave = zeros(length(tempo),1);
PMnave = zeros(length(tempo),1);
PMwave = zeros(length(tempo),1);

for i = 1:1:length(tempo)
    Cmave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*Cm(i,:)); %equation 20 supp applied to all for average calculations
    C2ave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*C2(i,:));
    C3ave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*C3(i,:));
    C4ave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*C4(i,:));
    C5ave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*C5(i,:));
    C6ave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*C6(i,:));
    C7ave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*C7(i,:));
    C8ave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*C8(i,:));
    C9ave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*C9(i,:));
    Cwave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*Cw(i,:));
    Cdsave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*Cds(i,:));
    Cdstot(i) = trapz(R, 4.*pi.*R.*R.*Cds(i,:));
    mu0ave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*mu0(i,:));
    mu1ave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*mu1(i,:));
    Mon(i) = trapz(R, 4.*pi.*R.*R.*mu1(i,:));
    mu2ave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*mu2(i,:));
    PMnave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*PMnum(i,:));
    PMwave(i) = 1/volume*trapz(R, 4.*pi.*R.*R.*PMw(i,:));
end

% Calculation of mass fluxes to evaluate mass loss
Kc_mon = Sh*Dmon_wat/(2*Radius);  %equation 10
Mass_loss = Mon./Mon(1);

Dmon_eff = Dmon_pol.*exp(beta.*sqrt(abs(1 - mu1./mu0.*PMmon./PMn))); %equation 15 paper

% Only for monomer/polymer
Cm_sup = Cm(:,Radius_Fragments)./(1 + Kc_mon*deltaR./Dmon_eff(:,Radius_Fragments));%equation 17 aproximation
C2_sup = C2(:,Radius_Fragments)./(1 + Kc_mon*deltaR./Dmon_eff(:,Radius_Fragments));
C3_sup = C3(:,Radius_Fragments)./(1 + Kc_mon*deltaR./Dmon_eff(:,Radius_Fragments));
C4_sup = C4(:,Radius_Fragments)./(1 + Kc_mon*deltaR./Dmon_eff(:,Radius_Fragments));
C5_sup = C5(:,Radius_Fragments)./(1 + Kc_mon*deltaR./Dmon_eff(:,Radius_Fragments));
C6_sup = C6(:,Radius_Fragments)./(1 + Kc_mon*deltaR./Dmon_eff(:,Radius_Fragments));
C7_sup = C7(:,Radius_Fragments)./(1 + Kc_mon*deltaR./Dmon_eff(:,Radius_Fragments));
C8_sup = C8(:,Radius_Fragments)./(1 + Kc_mon*deltaR./Dmon_eff(:,Radius_Fragments));
C9_sup = C9(:,Radius_Fragments)./(1 + Kc_mon*deltaR./Dmon_eff(:,Radius_Fragments));

Flux_mon = -Dmon_eff(:,Radius_Fragments).*(Cm_sup - Cm(:,Radius_Fragments))./deltaR.*SurfaceArea.*PMmon;
Flux_C2 = -Dmon_eff(:,Radius_Fragments).*(C2_sup - C2(:,Radius_Fragments))./deltaR.*SurfaceArea.*2.*PMmon;
Flux_C3 = -Dmon_eff(:,Radius_Fragments).*(C3_sup - C3(:,Radius_Fragments))./deltaR.*SurfaceArea.*3.*PMmon;
Flux_C4 = -Dmon_eff(:,Radius_Fragments).*(C4_sup - C4(:,Radius_Fragments))./deltaR.*SurfaceArea.*4.*PMmon;
Flux_C5 = -Dmon_eff(:,Radius_Fragments).*(C5_sup - C5(:,Radius_Fragments))./deltaR.*SurfaceArea.*5.*PMmon;
Flux_C6 = -Dmon_eff(:,Radius_Fragments).*(C6_sup - C6(:,Radius_Fragments))./deltaR.*SurfaceArea.*6.*PMmon;
Flux_C7 = -Dmon_eff(:,Radius_Fragments).*(C7_sup - C7(:,Radius_Fragments))./deltaR.*SurfaceArea.*7.*PMmon;
Flux_C8 = -Dmon_eff(:,Radius_Fragments).*(C8_sup - C8(:,Radius_Fragments))./deltaR.*SurfaceArea.*8.*PMmon;
Flux_C9 = -Dmon_eff(:,Radius_Fragments).*(C9_sup - C9(:,Radius_Fragments))./deltaR.*SurfaceArea.*9.*PMmon;
Flux_tot = Flux_mon + Flux_C2 + Flux_C3 + Flux_C4 + Flux_C5 + Flux_C6 + Flux_C7 + Flux_C8 + Flux_C9;
Mass_flux = cumtrapz(tempo, Flux_tot);
Mass_loss2 = (Massa - Mass_flux)/Massa;

Rel_sim = 1 - (Cdstot)./Cdstot(1);

% figure(1)
% plot(time_exp/86400, Rel_exp*100, 'o', tempo/86400, Rel_sim*100)
% title('Percentage Drug Release Over Time')
% xlabel('Time [days]')
% ylabel('Released drug [%]')
% legend('Experimental', 'Model')

% Print results in Excel file
h1 = {'Tempo [days]'};
h2 = {'Released drug [%]'};
h3 = {'Ddrug [cm2 s-1]'};

tempo = tempo/86400;
time_exp = time_exp/86400;

tempo = array2table(tempo, 'VariableNames', h1);
Rel_sim = array2table(Rel_sim, 'VariableNames', h2);
time_exp = array2table(time_exp, 'VariableNames', h1);
Rel_exp = array2table(Rel_exp, 'VariableNames', h2);
Ddrug_pol = array2table(Ddrug_pol, 'VariableNames', h2);

writetable(tempo,'Siepmann_et_al_sim_drug_release.xls','Range', 'A1')
writetable(Rel_sim,'Siepmann_et_al_sim_drug_release.xls','Range', 'B1')
writetable(time_exp,'Siepmann_et_al_sim_drug_release.xls','Range', 'C1')
writetable(Rel_exp,'Siepmann_et_al_sim_drug_release.xls','Range', 'D1')
writetable(Ddrug_pol,'Siepmann_et_al_sim_drug_release.xls','Range', 'E1')

tempo = table2array(tempo);
Rel_sim = table2array(Rel_sim);
time_exp = table2array(time_exp);
Rel_exp = table2array(Rel_exp);
Ddrug_pol = table2array(Ddrug_pol);


% writecell(h1,'Siepmann_et_al_sim_drug_release.xls', 'Range', 'A1')
% writecell(h2,'Siepmann_et_al_sim_drug_release.xls', 'Range', 'B1')
% writecell(h1,'Siepmann_et_al_sim_drug_release.xls', 'Range', 'C1')
% writecell(h2,'Siepmann_et_al_sim_drug_release.xls', 'Range', 'D1')
% writecell(h3,'Siepmann_et_al_sim_drug_release.xls', 'Range', 'E1')


% xlswrite('Siepmann_et_al_sim_drug_release.xls', tempo/86400, '106', 'A2');
% xlswrite('Siepmann_et_al_sim_drug_release.xls', Rel_sim, '106', 'B2');
% xlswrite('Siepmann_et_al_sim_drug_release.xls', time_exp/86400, '106', 'C2');
% xlswrite('Siepmann_et_al_sim_drug_release.xls', Rel_exp, '106', 'D2');
% xlswrite('Siepmann_et_al_sim_drug_release.xls', Ddrug_pol, '106', 'E2');



end

