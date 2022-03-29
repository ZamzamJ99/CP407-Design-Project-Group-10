close all
clear all

%Create graph to show variation of reaction rate with temperature

% Define general variables

M = 165.678675; % Ratio of inital concentrations of hydrogen (B) and oxygen (A), (CB0/CA0), (unitless);
s = 2 % stoichiometric ratio (b/a) (unitless)
CA0 = 0.1417445010928; %inital concentration of oxygen (A)(mol/m^3)
FA0 = 0.190931667; %  inital molar flowrate of Oxygen(mol/s) 
v0 = FA0/CA0; % volumetric flowrate of reactor (m^3/s)
Ti = 35:5:85; %Define temperature interval (^oC)
T = Ti+273.15; % express temperature interval in kelvin 
k0 = 1.96e11; %  pre-exponential factor (/s)
Ea = 73.77e3; % Activation energy (J/mol)
R = 8.3145; % Molar gas constant (J/mol K)
k = k0*exp((-Ea)./(R*T)); % rate constant (/s)
rA = k*CA0^2*((1-0.999996)*(M - (0.999996*s))); %reaction rate at varying temperatures (mol/m^3s)

% %Volume design
% % Define conversion interval
dXA = 0.999996/100;

% Define inital conditions 
XA(1) = 0; % set inital conversion of oxygen (A) as 0 (unitless)
V(1) = 0.05770956685; % set inital volume (m^3)

%Define variables
TR =  341.15; % Inlet operating temperature (K)
kR = k0*exp((-Ea)/(R*TR)); % operating rate constant (/s)

%iterate from XA = 0 to 0.999996 and calculate volume using Euler's method (V) and
%analytical solution (A)

for i = 2:101
    XA(i) = (i-1)*dXA;
    rAR(i) = kR.*CA0.^2.*(1 - XA(i-1)).*(M - (s*XA(i-1)));
    dVdXA = FA0./rAR(i);
    V(i) = V(i-1) + dVdXA*dXA;
    A(i) = (v0/(kR*CA0*(s-M)))*(log((1-XA(i))/(M-(s*XA(i))))-log(1/M))
end

% Find volume as predicted by both methods
V_Eu_R = V(:,end);
V_A_R = A(:,end);

%Plot reaction rate vs temperature
figure(1)
plot(Ti,rA)
xlabel('Temperature (^o C)')
ylabel('Reaction rate,-r_A,(mol / m^3s)')

%plot volume vs conversion
figure(2)
plot(XA,V,'r*-')
hold on
plot(XA,A,'b','LineWidth',1.1)
hold off
xlabel('Conversion of Oxygen (A), X_A')
ylabel('Volume, V, (m^3)')
legend('Euler method','Analytical Method')