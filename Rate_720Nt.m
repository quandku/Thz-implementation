clc;
clear all;
close all; 
%NantTX = [128,256,512];
%NantTX = [320,384,640];
NantTX = 720;
%% system parameters
Sysparam = struct();

%Sysparam.Nt = 256;
%Sysparam.N = Sysparam.Nt/Sysparam.M; 

Sysparam.Nr = 4;
Sysparam.Ns = 4;
Sysparam.fc = 3e11;%carrier frequency 
Sysparam.f = 30e9; %bandwidth 
Sysparam.K = 129;  %number of sub-carriers
Sysparam.M = 16; %number of TTDs for each RF chain
Sysparam.L = 4; % number of paths equals to number of RF chains
% compute discrete time-delay values
Sysparam.Q = 16; % number of fixed TTD element 
Sysparam.tmax = (1e-12)*1000;
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
Simparam.phi = 0.2951*pi*ones(Sysparam.L,1);%AoD;
%Simparam.phi = pi*rand(Sysparam.L,1)-0.5*pi;%AoD;
%% Diferent transmit power
%P = -5:5:5;% transmit power in dB
P = 0;
Rho = 10.^(P/10); %transmit power
%P = 1; 

%% Achievable rate comparison
Rate1 = zeros(Sysparam.K,length(NantTX)); % ideal Analog precoder
Rate2 = zeros(Sysparam.K,length(NantTX)); % joint TTD and PS 
Rate3 = zeros(Sysparam.K,length(NantTX)); % fixed PS 
Rate4 = zeros(Sysparam.K,length(NantTX)); % Dynamic subarray
Nrlz = 100; %number of channel realizations
Nloop = 30;
%% 
for n =1:length(NantTX)
Sysparam.Nt = NantTX(n);
Sysparam.N = Sysparam.Nt/Sysparam.M; 

%% Ideal Analog precoder
for l = 1:Sysparam.L
    %F3(:,l) = ArrResponseGenerateH(Sysparam.Nt,Simparam.phi(l),1)/sqrt(Sysparam.Nt);
    for k =1:Sysparam.K
        idealF(:,l,k) = ArrResponseGenerateH(Sysparam.Nt,Simparam.phi(l),Sysparam.xi(k))/sqrt(Sysparam.Nt);
    end 
end 
for nRlz = 1:Nrlz
    %% generate channel
    H = WidebandChannel(Sysparam,Simparam);   
    %% compensate beam squint
    [G1,G2k,~,~] = fixPSdesignTTD(Sysparam,Simparam); %Tan design
    [PS,TTD,~,~] = jointPSandTTD(Sysparam, Simparam); % joint delay and phase precoding
   [S,F2k,D,~] = DynamicSubarray(Sysparam,Simparam,H,Nloop); %Dynamic subarray method
    %%  
    for k = 1:Sysparam.K
       %% ideal analog precoder
          Heff1(:,:,k) = H(:,:,k)*idealF(:,:,k);
          [~,~,Veff1(:,:,k)] = svd(Heff1(:,:,k));
          Rate1(k,n) = Rate1(k,n)+log2(det(eye(Sysparam.Ns)+(Rho/Sysparam.Ns)*Heff1(:,:,k)*Veff1(:,1:Sysparam.Ns,k)*Veff1(:,1:Sysparam.Ns,k)'*Heff1(:,:,k)'));         
       %% Joint PS and TTD rate - Rate 2   
          Heff2(:,:,k) = H(:,:,k)*PS*TTD(:,:,k);
          [~,~,Veff2(:,:,k)] = svd(Heff2(:,:,k));
          Rate2(k,n) = Rate2(k,n)+log2(det(eye(Sysparam.Ns)+(Rho/Sysparam.Ns)*Heff2(:,:,k)*Veff2(:,1:Sysparam.Ns,k)*Veff2(:,1:Sysparam.Ns,k)'*Heff2(:,:,k)'));         
         %% Fixed PS rate - Rate 3   
         Heff3(:,:,k) = H(:,:,k)*G1*G2k(:,:,k);
         [~,~,Veff3(:,:,k)] = svd(Heff3(:,:,k));
         Rate3(k,n) = Rate3(k,n)+log2(det(eye(Sysparam.Ns)+(Rho/Sysparam.Ns)*Heff3(:,:,k)*Veff3(:,1:Sysparam.Ns,k)*Veff3(:,1:Sysparam.Ns,k)'*Heff3(:,:,k)'));
         %% Dynamic subarray rate - Rate 4
         Heff4(:,:,k) = H(:,:,k)*S*F2k(:,:,k);
         [~,~,Veff4(:,:,k)] = svd(Heff4(:,:,k));
         Rate4(k,n) = Rate4(n)+log2(det(eye(Sysparam.Ns)+(Rho/Sysparam.Ns)*Heff4(:,:,k)*Veff4(:,1:Sysparam.Ns,k)*Veff4(:,1:Sysparam.Ns,k)'*Heff4(:,:,k)'));    
    end  
end
clear idealF 
end 
%%
Rate1 = real(Rate1)./Nrlz;
Rate2 = real(Rate2)./Nrlz;
Rate3 = real(Rate3)./Nrlz;
Rate4 = real(Rate4)./Nrlz;
%% 
% 
% AverageRate1 = sum(real(Rate1))./(Nrlz*Sysparam.K)
% AverageRate2 = sum(real(Rate2))./(Nrlz*Sysparam.K)
% AverageRate3 = sum(real(Rate3))./(Nrlz*Sysparam.K)
%%
set(0, 'defaultlinelinewidth', 1); set(0, 'defaultlinemarkersize', 8);
set(0, 'defaultaxesfontsize', 15); set(0, 'defaulttextfontsize', 15); 
for n = 1:length(NantTX)
figure(n)
f2 = cdfplot(Rate2(:,n));
set(f2,'LineWidth',4,'LineStyle',':','Marker','none');
hold on 
f3 = cdfplot(Rate3(:,n));
set(f3,'LineWidth',4,'LineStyle','--','Marker','none');
% f4 = cdfplot(Rate4(:,n));
% set(f4,'LineWidth',4,'LineStyle','--','Marker','none');
f1 = cdfplot(Rate1(:,n));
set(f1,'LineWidth',4,'LineStyle','-','Marker','none');
hold on
legend('Proposed approach','TTD-based precoding [8]','Ideal analog precoder')
xlabel('Per-subcarrier achievable rate (bit/s/Hz)')
ylabel('CDF')
title('')
end 
