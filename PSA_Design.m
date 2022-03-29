close all
clear all

P = 40e3; %Operating pressure (mBar)
T = (30+273.15); % Temperature (K)
e = 0.4; % bed void fraction (unitless)
u = 0.2679; % intersitial velocity (m/s)
a = 0.0375; % Kinetic driving factor (/s)
qm = 0.285/0.01802; % saturated loading (mol/kg)
yi = 0.00020385186910942; % mol fraction in inlet (unitless)
p0 = yi*P; %partial pressure (mbar)
rhop = 1000; % water density (kg/m3)
b = 2.03629e-3; %(/mbar)
R = 8.3145; %Ideal gas constant 

% Define length interval
L = 5; %set upper bound for length
N = 101; % set number of iterations/distance between length
z = linspace(0,L,N); %create length variable
dz = z(2)-z(1); %calculate change in length

% Define time interval
t = 0:1800

%inital condition
ICA = zeros(1,N);
ICB = zeros(1,N);
IC = [ICA,ICB]; %create variable for inital conditions of pressure (p) and loading (q)

%ODE solver
[t y] = ode15s(@f,t,IC,[],e,u,a,b,qm,p0,rhop,N,dz,R,T); %solve function and store it in these variables

%extract result
p = y(:,1:N); 
q = y(:,N+1:2*N);

%reinput boundary condition
p(:,1) = p0;
p(:,end) = (4*p(:,end-1)- p(:,end-2))./3; %use finite difference method to compute p and q using final condition
q(:,end) = (4*q(:,end-1)- q(:,end-2))./3;

figure(1)
imagesc(z,t,p)
colorbar
grid on
ylabel('Time (s)')
xlabel('Height of bed (m)')

p_tc = p(600,:); % desired cycle time tc = 600s
p_idx = find(p_tc >0.15 & p_tc<0.161,1); %find 1st value of desired purity at cycle time
pp = p(600,p_idx); %check index of 1st value is desired purity
z_op = z(p_idx); %Find bed operating height at desired cycle time and purity
p_zop = p(:,p_idx); %find pressure profile at this height

%conduct analysis on how height will affect composition
p_zidx = p_idx + 2; % increase height by 0.1m
p_zi = p(:,p_zidx); %find pressure profile at this height
p_tidx = find(p_zi>0.15 & p_zi<0.161,1); %find time index of 1st value of desired purity
pp2 = p(p_tidx, p_zidx); %check value of purity is the same as earlier

tc_ip = ((p_tidx/600)-1)*100; % find percentage cycle time increase per 0.1m height increase

p_2x = ((0.1512/0.1502)-1) % increase when composition is doubled
n_i = ((0.16/0.1502) - 1)/p_2x % number of times composition can be doubled
yi_max_ppm = n_i*2*yi*1e6; % find maximum inlet mole fraction of water vapour PSA can treat

figure(2)
%plot breakthrough curve at operating height
pplot = p_zop./p0;
plot(t,pplot,'b')
xline(600,'r')
xlabel('Time (s)')
ylabel('p /p_0 (unitless)')


%define function to solve
function dydt = f(t,y,e,u,a,b,qm,p0,rhop,N,dz,R,T)
dydt = zeros(length(y),1);
dpdt = zeros(N,1);
dqdt = zeros(N,1);

% define value

p = y(1:N);
q = y(N+1:2*N);

%boundary condition
p(1) = p0;
p(end) = (4*p(end-1)- p(end-2))./3 % use finite difference method to estimate p and q using final condition

for i = 2:N-1 
    dqdt(i) = a.*((qm.*b.*p(i)./(1+(b.*p(i))))-q(i));
    dpdz(i) = (p(i+1)-p(i-1))./(2*dz);
    dpdt(i) = (-u.*(dpdz(i))) - (R.*T.*((1-e)./e).*dqdt(i));
end

dydt = [dpdt;dqdt]; 
end
    
