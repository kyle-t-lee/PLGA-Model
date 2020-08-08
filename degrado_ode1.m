function dy = degrado_ode1(t,y)

global Radius Radius_Fragments Kd deltaR Sh Dmon_pol Dmon_wat Dwat_wat Dwat_pol Cw_bulk
global PMmon PMn beta Ddrug_wat Ddrug_pol

Neq = 14*Radius_Fragments;

dy = zeros(Neq,1);

R = deltaR:deltaR:Radius;

% Variables initialization
for i = 1:1:Radius_Fragments
   Cm(i) = y(i);
   C2(i) = y(i + Radius_Fragments);
   C3(i) = y(i + 2*Radius_Fragments);
   C4(i) = y(i + 3*Radius_Fragments);
   C5(i) = y(i + 4*Radius_Fragments);
   C6(i) = y(i + 5*Radius_Fragments);
   C7(i) = y(i + 6*Radius_Fragments);
   C8(i) = y(i + 7*Radius_Fragments);
   C9(i) = y(i + 8*Radius_Fragments);
   Cw(i) = y(i + 9*Radius_Fragments);
   mu0(i) = y(i + 10*Radius_Fragments);
   mu1(i) = y(i + 11*Radius_Fragments);
   mu2(i) = y(i + 12*Radius_Fragments);
   Cds(i) = y(i + 13*Radius_Fragments);
end

%equation 15 paper
PMnum = mu1./mu0.*PMmon;
dPMn = abs(1 - PMnum./PMn);

Dmon_eff = Dmon_pol.*exp(beta.*sqrt(dPMn));
Dwat_eff = Dwat_pol.*exp(beta.*sqrt(dPMn));
Ddrug_eff = Ddrug_pol.*exp(beta.*sqrt(dPMn));


% Boundary conditions equation 10
Kc_mon = Sh*Dmon_wat/(2*Radius);
Kc_wat = Sh*Dwat_wat/(2*Radius);
Kc_drug = Sh*Ddrug_wat/(2*Radius);

Cm_sup = Cm(Radius_Fragments)/(1 + Kc_mon*deltaR/Dmon_eff(Radius_Fragments));
C2_sup = C2(Radius_Fragments)/(1 + Kc_mon*deltaR/Dmon_eff(Radius_Fragments));
C3_sup = C3(Radius_Fragments)/(1 + Kc_mon*deltaR/Dmon_eff(Radius_Fragments));
C4_sup = C4(Radius_Fragments)/(1 + Kc_mon*deltaR/Dmon_eff(Radius_Fragments));
C5_sup = C5(Radius_Fragments)/(1 + Kc_mon*deltaR/Dmon_eff(Radius_Fragments));
C6_sup = C6(Radius_Fragments)/(1 + Kc_mon*deltaR/Dmon_eff(Radius_Fragments));
C7_sup = C7(Radius_Fragments)/(1 + Kc_mon*deltaR/Dmon_eff(Radius_Fragments));
C8_sup = C8(Radius_Fragments)/(1 + Kc_mon*deltaR/Dmon_eff(Radius_Fragments));
C9_sup = C9(Radius_Fragments)/(1 + Kc_mon*deltaR/Dmon_eff(Radius_Fragments));
Cw_sup = (Cw(Radius_Fragments) + Cw_bulk*Kc_wat*deltaR/Dwat_eff(Radius_Fragments))/(1 + Kc_wat*deltaR/Dwat_eff(Radius_Fragments));
Cdrug_sup = Cds(Radius_Fragments)/(1 + Kc_drug*deltaR/Ddrug_eff(Radius_Fragments));

for i = 1:1:Radius_Fragments
    if i == 1
        der1Cm(i) = 0; %Equation 8 supplimetary
        der1C2(i) = 0; %Equation 8 supplimetary
        der1C3(i) = 0; %Equation 8 supplimetary
        der1C4(i) = 0; %Equation 8 supplimetary
        der1C5(i) = 0; %Equation 8 supplimetary
        der1C6(i) = 0; %Equation 8 supplimetary
        der1C7(i) = 0; %Equation 8 supplimetary
        der1C8(i) = 0; %Equation 8 supplimetary
        der1C9(i) = 0; %Equation 8 supplimetary
        der1Cw(i) = 0; %Equation 8 supplimetary
        der1Dm(i) = 0; %Equation 8 supplimetary
        der1Dw(i) = 0; %Equation 8 supplimetary
        der1Ddrug(i) = 0;  %Equation 8 supplimetary
        der1Cds(i) = 0; %Equation 8 supplimetary
        der2Cm(i) = (Cm(i)-2*Cm(i)+Cm(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C2(i) = (C2(i)-2*C2(i)+C2(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C3(i) = (C3(i)-2*C3(i)+C3(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C4(i) = (C4(i)-2*C4(i)+C4(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C5(i) = (C5(i)-2*C5(i)+C5(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C6(i) = (C6(i)-2*C6(i)+C6(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C7(i) = (C7(i)-2*C7(i)+C7(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C8(i) = (C8(i)-2*C8(i)+C8(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C9(i) = (C9(i)-2*C9(i)+C9(i+1))/(deltaR^2); % numerical expression for second derivative
        der2Cw(i) = (Cw(i)-2*Cw(i)+Cw(i+1))/(deltaR^2); % numerical expression for second derivative
        der2Cds(i) = (Cds(i)-2*Cds(i)+Cds(i+1))/(deltaR^2); % numerical expression for second derivative
    elseif i == Radius_Fragments
        der1Cm(i) = (Cm_sup-Cm(i))/deltaR;  %Equation 9 supplimetary
        der1C2(i) = (C2_sup-C2(i))/deltaR;  %Equation 9 supplimetary
        der1C3(i) = (C3_sup-C3(i))/deltaR;  %Equation 9 supplimetary
        der1C4(i) = (C4_sup-C4(i))/deltaR;  %Equation 9 supplimetary
        der1C5(i) = (C5_sup-C5(i))/deltaR;  %Equation 9 supplimetary
        der1C6(i) = (C6_sup-C6(i))/deltaR;  %Equation 9 supplimetary
        der1C7(i) = (C7_sup-C7(i))/deltaR;  %Equation 9 supplimetary
        der1C8(i) = (C8_sup-C8(i))/deltaR;  %Equation 9 supplimetary
        der1C9(i) = (C9_sup-C9(i))/deltaR;  %Equation 9 supplimetary
        der1Cw(i) = (Cw_sup-Cw(i))/deltaR;  %Equation 9 supplimetary
        der1Cds(i) = (Cdrug_sup-Cds(i))/deltaR;  %Equation 9 supplimetary
        der1Dm(i) = 0;  %Equation 9 supplimetary
        der1Dw(i) = 0;  %Equation 9 supplimetary
        der1Ddrug(i) = 0;  %Equation 9 supplimetary
        der2Cm(i) = (Cm(i-1)-2*Cm(i)+Cm_sup)/(deltaR^2); % numerical expression for second derivative
        der2C2(i) = (C2(i-1)-2*C2(i)+C2_sup)/(deltaR^2); % numerical expression for second derivative
        der2C3(i) = (C3(i-1)-2*C3(i)+C3_sup)/(deltaR^2); % numerical expression for second derivative
        der2C4(i) = (C4(i-1)-2*C4(i)+C4_sup)/(deltaR^2); % numerical expression for second derivative
        der2C5(i) = (C5(i-1)-2*C5(i)+C5_sup)/(deltaR^2); % numerical expression for second derivative
        der2C6(i) = (C6(i-1)-2*C6(i)+C6_sup)/(deltaR^2); % numerical expression for second derivative
        der2C7(i) = (C7(i-1)-2*C7(i)+C7_sup)/(deltaR^2); % numerical expression for second derivative
        der2C8(i) = (C8(i-1)-2*C8(i)+C8_sup)/(deltaR^2); % numerical expression for second derivative
        der2C9(i) = (C9(i-1)-2*C9(i)+C9_sup)/(deltaR^2); % numerical expression for second derivative
        der2Cw(i) = (Cw(i-1)-2*Cw(i)+Cw_sup)/(deltaR^2); % numerical expression for second derivative
        der2Cds(i) = (Cds(i-1)-2*Cds(i)+Cdrug_sup)/(deltaR^2); % numerical expression for second derivative
    else
        der1Cm(i) = (Cm(i+1)-Cm(i-1))/(2*deltaR); % numerical expression for first derivative
        der1C2(i) = (C2(i+1)-C2(i-1))/(2*deltaR); % numerical expression for first derivative
        der1C3(i) = (C3(i+1)-C3(i-1))/(2*deltaR); % numerical expression for first derivative
        der1C4(i) = (C4(i+1)-C4(i-1))/(2*deltaR); % numerical expression for first derivative
        der1C5(i) = (C5(i+1)-C5(i-1))/(2*deltaR); % numerical expression for first derivative
        der1C6(i) = (C6(i+1)-C6(i-1))/(2*deltaR); % numerical expression for first derivative
        der1C7(i) = (C7(i+1)-C7(i-1))/(2*deltaR); % numerical expression for first derivative
        der1C8(i) = (C8(i+1)-C8(i-1))/(2*deltaR); % numerical expression for first derivative
        der1C9(i) = (C9(i+1)-C9(i-1))/(2*deltaR); % numerical expression for first derivative
        der1Cds(i) = (Cds(i+1)-Cds(i-1))/(2*deltaR); % numerical expression for first derivative
        der1Cw(i) = (Cw(i+1)-Cw(i-1))/(2*deltaR); % numerical expression for first derivative
        der1Dm(i) = (Dmon_eff(i+1)-Dmon_eff(i-1))/(2*deltaR); % numerical expression for first derivative
        der1Dw(i) = (Dwat_eff(i+1)-Dwat_eff(i-1))/(2*deltaR); % numerical expression for first derivative
        der1Ddrug(i) = (Ddrug_eff(i+1)-Ddrug_eff(i-1))/(2*deltaR); % numerical expression for first derivative
        der2Cm(i) = (Cm(i-1)-2*Cm(i)+Cm(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C2(i) = (C2(i-1)-2*C2(i)+C2(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C3(i) = (C3(i-1)-2*C3(i)+C3(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C4(i) = (C4(i-1)-2*C4(i)+C4(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C5(i) = (C5(i-1)-2*C5(i)+C5(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C6(i) = (C6(i-1)-2*C6(i)+C6(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C7(i) = (C7(i-1)-2*C7(i)+C7(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C8(i) = (C8(i-1)-2*C8(i)+C8(i+1))/(deltaR^2); % numerical expression for second derivative
        der2C9(i) = (C9(i-1)-2*C9(i)+C9(i+1))/(deltaR^2); % numerical expression for second derivative
        der2Cds(i) = (Cds(i-1)-2*Cds(i)+Cds(i+1))/(deltaR^2); % numerical expression for second derivative
        der2Cw(i) = (Cw(i-1)-2*Cw(i)+Cw(i+1))/(deltaR^2); % numerical expression for second derivative
    end
    
    DiffC1(i) = Dmon_eff(i)*(2/R(i)*der1Cm(i) + der2Cm(i)) + der1Dm(i)*der1Cm(i); % equation 7 supplimetary
    DiffC2(i) = Dmon_eff(i)*(2/R(i)*der1C2(i) + der2C2(i)) + der1Dm(i)*der1C2(i); % equation 7 supplimetary
    DiffC3(i) = Dmon_eff(i)*(2/R(i)*der1C3(i) + der2C3(i)) + der1Dm(i)*der1C3(i); % equation 7 supplimetary
    DiffC4(i) = Dmon_eff(i)*(2/R(i)*der1C4(i) + der2C4(i)) + der1Dm(i)*der1C4(i); % equation 7 supplimetary
    DiffC5(i) = Dmon_eff(i)*(2/R(i)*der1C5(i) + der2C5(i)) + der1Dm(i)*der1C5(i); % equation 7 supplimetary
    DiffC6(i) = Dmon_eff(i)*(2/R(i)*der1C6(i) + der2C6(i)) + der1Dm(i)*der1C6(i); % equation 7 supplimetary
    DiffC7(i) = Dmon_eff(i)*(2/R(i)*der1C7(i) + der2C7(i)) + der1Dm(i)*der1C7(i); % equation 7 supplimetary
    DiffC8(i) = Dmon_eff(i)*(2/R(i)*der1C8(i) + der2C8(i)) + der1Dm(i)*der1C8(i); % equation 7 supplimetary
    DiffC9(i) = Dmon_eff(i)*(2/R(i)*der1C9(i) + der2C9(i)) + der1Dm(i)*der1C9(i); % equation 7 supplimetary
    DiffCds(i) = Ddrug_eff(i)*(2/R(i)*der1Cds(i) + der2Cds(i)) + der1Ddrug(i)*der1Cds(i); % equation 7 supplimetary
        
    dy(i) = DiffC1(i) + 2*Kd*Cw(i)*(mu0(i)-Cm(i))*mu0(i); %Equation 5.a
    dy(i + Radius_Fragments) = DiffC2(i) + 2*Kd*Cw(i)*(mu0(i)-Cm(i)-C2(i))*mu0(i) - (2-1)*Kd*Cw(i)*C2(i)*mu0(i);%Equation 5.b n=2
    dy(i + 2*Radius_Fragments) = DiffC3(i) + 2*Kd*Cw(i)*(mu0(i)-Cm(i)-C2(i)-C3(i))*mu0(i) - (3-1)*Kd*Cw(i)*C3(i)*mu0(i);%Equation 5.b n=3
    dy(i + 3*Radius_Fragments) = DiffC4(i) + 2*Kd*Cw(i)*(mu0(i)-Cm(i)-C2(i)-C3(i)-C4(i))*mu0(i) - (4-1)*Kd*Cw(i)*C4(i)*mu0(i);%Equation 5.b n=4
    dy(i + 4*Radius_Fragments) = DiffC5(i) + 2*Kd*Cw(i)*(mu0(i)-Cm(i)-C2(i)-C3(i)-C4(i)-C5(i))*mu0(i) - (5-1)*Kd*Cw(i)*C5(i)*mu0(i);%Equation 5.b n5
    dy(i + 5*Radius_Fragments) = DiffC6(i) + 2*Kd*Cw(i)*(mu0(i)-Cm(i)-C2(i)-C3(i)-C4(i)-C5(i)-C6(i))*mu0(i) - (6-1)*Kd*Cw(i)*C6(i)*mu0(i);%Equation 5.b n=6
    dy(i + 6*Radius_Fragments) = DiffC7(i) + 2*Kd*Cw(i)*(mu0(i)-Cm(i)-C2(i)-C3(i)-C4(i)-C5(i)-C6(i)-C7(i))*mu0(i) - (7-1)*Kd*Cw(i)*C7(i)*mu0(i);%Equation 5.b n=7
    dy(i + 7*Radius_Fragments) = DiffC8(i) + 2*Kd*Cw(i)*(mu0(i)-Cm(i)-C2(i)-C3(i)-C4(i)-C5(i)-C6(i)-C7(i)-C8(i))*mu0(i) - (8-1)*Kd*Cw(i)*C8(i)*mu0(i);%Equation 5.b n=8
    dy(i + 8*Radius_Fragments) = DiffC9(i) + 2*Kd*Cw(i)*(mu0(i)-Cm(i)-C2(i)-C3(i)-C4(i)-C5(i)-C6(i)-C7(i)-C8(i)-C9(i))*mu0(i) - (9-1)*Kd*Cw(i)*C9(i)*mu0(i);%Equation 5.b n=9
    dy(i + 9*Radius_Fragments) = Dwat_eff(i)*(2/R(i)*der1Cw(i) + der2Cw(i)) + der1Dw(i)*der1Cw(i) - Kd*Cw(i)*(mu1(i)-mu0(i))*mu0(i);%Equation 5.c
    dy(i + 10*Radius_Fragments) = DiffC1(i) + DiffC2(i) + DiffC3(i) + DiffC4(i) + DiffC5(i) + DiffC6(i) + DiffC7(i) + DiffC8(i) + DiffC9(i) + Kd*Cw(i)*(mu1(i)-mu0(i))*mu0(i);%Equation 5.d
    dy(i + 11*Radius_Fragments) = DiffC1(i) + 2*DiffC2(i) + 3*DiffC3(i) + 4*DiffC4(i) + 5*DiffC5(i) + 6*DiffC6(i) + 7*DiffC7(i) + 8*DiffC8(i) + 9*DiffC9(i);%Equation 5.e
    dy(i + 12*Radius_Fragments) = DiffC1(i) + 4*DiffC2(i) + 9*DiffC3(i) + 16*DiffC4(i) + 25*DiffC5(i) + 36*DiffC6(i) + 49*DiffC7(i) + 64*DiffC8(i) + 81*DiffC9(i) + Kd/3*Cw(i)*(mu1(i)-2*(mu2(i)^2)/mu1(i) + mu2(i)*mu1(i)/mu0(i))*mu0(i);%Equation 5.f    
    dy(i + 13*Radius_Fragments) = DiffCds(i);
end