clc;
clear all;
close all; 
%% system parameters
Sysparam = struct();

Sysparam.Nt = 256;
Sysparam.Nr = 4;
Sysparam.Ns = 4;
Sysparam.fc = 3e11;%carrier frequency 
Sysparam.f = 30e9; %bandwidth 
Sysparam.K = 129;  %number of sub-carriers
Sysparam.M = 8; %number of TTDs for each RF chain
Sysparam.N = Sysparam.Nt/Sysparam.M;
Sysparam.L = 1 ; % number of paths equals to number of RF chains
% compute discrete time-delay values
Sysparam.Q = 8; % number of fixed TTD element 
Sysparam.tmax = (1e-12)*340;
Sysparam.T = Sysparam.tmax*(0:1:Sysparam.Q-1)/(Sysparam.Q-1);

% Sub-carrier frequencies
for k =1:Sysparam.K 
    Sysparam.F(k,1) = Sysparam.fc +(Sysparam.f/Sysparam.K)*(k-1-(Sysparam.K-1)/2); % sub-carrier frequency
end 

% relative frequencies
Sysparam.xi = Sysparam.F/Sysparam.fc; 

%% simulation parameters
Simparam = struct();
Simparam.psi = pi*rand(Sysparam.L,1)-0.5*pi; %AoA;
Simparam.phi = 0.2951*pi*ones(Sysparam.L,1)%AoD;
%% Diferent transmit power
P = -5:5:5;% transmit power in dB
Rho = 10.^(P/10); %transmit power
%P = 1; 
%% Analog precoder without beam squint compensation
for l = 1:Sysparam.L
    F3(:,l) = ArrResponseGenerateH(Sysparam.Nt,Simparam.phi(l),1)/sqrt(Sysparam.Nt);
end 
%% Achievable rate comparison
Rate1 = zeros(Sysparam.K,length(P)); %joint delay and phase precoding
Rate2 = zeros(Sysparam.K,length(P)); %fixed PS and design TTD 
Rate3 = zeros(Sysparam.K,length(P)); % fixed TTD and design switch network
Rate4 = zeros(Sysparam.K,length(P)); % ideal analog precoder
Rate5 = zeros(Sysparam.K,length(P)); % without beam squint compensation
% use to run BCD for dynamic subarray method
Nloop = 20;
Nrlz = 1; %number of channel realization
%%
%% 
for n = 1:length(P)
    for nRlz = 1:Nrlz
    %% generate channel
    H = WidebandChannel(Sysparam,Simparam);   
    %% compensate beam squint
    [G1,G2k,~,~] = fixPSdesignTTD(Sysparam,Simparam); %Tan design
    [PS,TTD,~,~] = jointPSandTTD(Sysparam, Simparam); % joint delay and phase precoding
    [S,F2k,D,~] = DynamicSubarray(Sysparam,Simparam,H,Nloop); %Dynamic subarray method
    %%  
    for k = 1:Sysparam.K
        %% Fully digital optimal rate -Rate 1
          [~,~,V(:,:,k)] = svd(H(:,:,k));
          Rate1(k,n) = Rate1(k,n)+log2(det(eye(Sysparam.Ns)+(Rho(n)/Sysparam.Ns)*H(:,:,k)*V(:,1:Sysparam.Ns,k)*V(:,1:Sysparam.Ns,k)'*H(:,:,k)'));
        %% Joint PS and TTD rate - Rate 2   
%           Heff2(:,:,k) = H(:,:,k)*PS*TTD(:,:,k);
%           [~,~,Veff2(:,:,k)] = svd(Heff2(:,:,k));
%           Rate2(k,n) = Rate2(k,n)+log2(det(eye(Sysparam.Ns)+(Rho(n)/Sysparam.Ns)*Heff2(:,:,k)*Veff2(:,1:Sysparam.Ns,k)*Veff2(:,1:Sysparam.Ns,k)'*Heff2(:,:,k)'));
%          %% Fixed PS and design TTD rate - Rate 3
%          Heff3(:,:,k) = H(:,:,k)*G1*G2k(:,:,k);
%          [~,~,Veff3(:,:,k)] = svd(Heff3(:,:,k));
%          Rate3(k,n) = Rate3(k,n)+log2(det(eye(Sysparam.Ns)+(Rho(n)/Sysparam.Ns)*Heff3(:,:,k)*Veff3(:,1:Sysparam.Ns,k)*Veff3(:,1:Sysparam.Ns,k)'*Heff3(:,:,k)'));
%          %% Dynamic subarray rate - Rate 4
%          Heff4(:,:,k) = H(:,:,k)*S*F2k(:,:,k);
%          [~,~,Veff4(:,:,k)] = svd(Heff4(:,:,k));
%          Rate4(k,n) = Rate4(k,n)+log2(det(eye(Sysparam.Ns)+(Rho(n)/Sysparam.Ns)*Heff4(:,:,k)*Veff4(:,1:Sysparam.Ns,k)*Veff4(:,1:Sysparam.Ns,k)'*Heff4(:,:,k)'));    
%          %% Without beam squint compensation
%          Heff5(:,:,k) = H(:,:,k)*F3;
%          [~,~,Veff5(:,:,k)] = svd(Heff5(:,:,k));
%          Rate5(k,n) = Rate5(k,n)+log2(det(eye(Sysparam.Ns)+(Rho(n)/Sysparam.Ns)*Heff5(:,:,k)*Veff5(:,1:Sysparam.Ns,k)*Veff5(:,1:Sysparam.Ns,k)'*Heff5(:,:,k)'));
    end  
    end 
end