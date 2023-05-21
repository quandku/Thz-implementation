close all;
clear variables;
%% system parameters
Sysparam = struct();

Sysparam.Nt = 256;
Sysparam.Nr = 4;
Sysparam.Nrf=4;
Sysparam.Ns = 4;
Sysparam.fc = 300*1e9;%carrier frequency 
Sysparam.f = 30*1e9; %bandwidth 
Sysparam.K = 129;  %number of sub-carriers
Sysparam.M = 4; %number of TTDs for each RF chain
Sysparam.N = Sysparam.Nt/Sysparam.M; %number of PSs for each TTD
Sysparam.L = 4 ; % number of paths equals to number of RF chains
Sysparam.tmax = 100e-12;% maximum time-delay value 
%% simulation parameters
Simparam = struct();

Simparam.phi = pi*rand(Sysparam.L,1)-0.5*pi; %AoD;
Simparam.psi = pi*rand(Sysparam.L,1)-0.5*pi; %AoA;

Simparam.P = -20:5:15;% transmit power in dB
Simparam.Rho = 10.^(Simparam.P/10); %transmit power
Simparam.Niter = 50; %channel realizations
%% Sub-carrier frequencies
for k =1:Sysparam.K 
    Sysparam.F(k,1) = Sysparam.fc +(Sysparam.f/Sysparam.K)*(k-1-(Sysparam.K-1)/2); % sub-carrier frequency
end 
% relative frequencies
Sysparam.xi = Sysparam.F/Sysparam.fc; 
%% Beam squint compensation performance
[PS,TTD,x,t] = jointPSandTTD(Sysparam,Simparam);
