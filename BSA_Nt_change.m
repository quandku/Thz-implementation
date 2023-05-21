clc;
clear all;
close all; 
NantTX = [128,256,512];
%NantTX = [320,384,640];
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
Sysparam.Q = 64; % number of fixed TTD element 
Sysparam.tmax = (1e-12)*340;
Sysparam.T = Sysparam.tmax*(0:1:Sysparam.Q-1)/(Sysparam.Q-1);
% Sub-carrier frequencies
for k =1:Sysparam.K 
    Sysparam.F(k,1) = Sysparam.fc +(Sysparam.f/Sysparam.K)*(k-1-(Sysparam.K-1)/2); % sub-carrier frequency
end 
% relative frequencies
Sysparam.xi = Sysparam.F/Sysparam.fc; 

%% simulation parameters
% Simparam = struct();
% Simparam.psi = pi*rand(Sysparam.L,1)-0.5*pi; %AoA;
% Simparam.phi = 0.2951*pi*ones(Sysparam.L,1);%AoD;
% SD = [0.65,0.75,0.8,0.85];
% Simparam.phi = asin(SD)./pi;
% Simparam.phi = pi*((1/6)*rand(Sysparam.L,1)+1/3);%AoD;
% % spatialDirection = Simparam.phi; 
% % sin(Simparam.phi)
load SimParamPhi.mat
%% Diferent transmit power
%P = -5:5:5;% transmit power in dB
P = 3;
Rho = 10.^(P/10); %transmit power
%P = 1; 

%% Achievable rate comparison
Rate1 = zeros(Sysparam.K,length(NantTX)); % ideal Analog precoder
Rate2 = zeros(Sysparam.K,length(NantTX)); % joint TTD and PS 
Rate3 = zeros(Sysparam.K,length(NantTX)); % fixed PS 
Rate4 = zeros(Sysparam.K,length(NantTX)); % Dynamic subarray
Rate5 = zeros(Sysparam.K,length(NantTX)); % Beam squint aware
Rate6 = zeros(Sysparam.K,length(NantTX)); % Without beam squint compensation
Nrlz = 100; %number of channel realizations
Nloop = 50;
%% 
for n =1:length(NantTX)
Sysparam.Nt = NantTX(n);
Sysparam.N = Sysparam.Nt/Sysparam.M; 

%% Ideal Analog precoder
for l = 1:Sysparam.L
    F3(:,l) = ArrResponseGenerateH(Sysparam.Nt,Simparam.phi(l),1)/sqrt(Sysparam.Nt);
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
    %%  
    for k = 1:Sysparam.K
       %% Joint PS and TTD rate - Rate 2   
        Heff2(:,:,k) = H(:,:,k)*PS*TTD(:,:,k);
        [~,~,Veff2(:,:,k)] = svd(Heff2(:,:,k));
        Rate2(k,n) = Rate2(k,n)+log2(det(eye(Sysparam.Ns)+(Rho/Sysparam.Ns)*Heff2(:,:,k)*Veff2(:,1:Sysparam.Ns,k)*Veff2(:,1:Sysparam.Ns,k)'*Heff2(:,:,k)'));         
        %% Fixed PS rate - Rate 3   
         Heff3(:,:,k) = H(:,:,k)*G1*G2k(:,:,k);
         [~,~,Veff3(:,:,k)] = svd(Heff3(:,:,k));
         Rate3(k,n) = Rate3(k,n)+log2(det(eye(Sysparam.Ns)+(Rho/Sysparam.Ns)*Heff3(:,:,k)*Veff3(:,1:Sysparam.Ns,k)*Veff3(:,1:Sysparam.Ns,k)'*Heff3(:,:,k)'));
        %% Beam squint aware
        Heff5(:,:,k) = H(:,:,k)*F3;
        [~,~,Veff5(:,:,k)] = svd(Heff5(:,:,k));
        V5(:,:,k) = pinv(F3)*idealF(:,:,k)*Veff5(:,:,k);%beam squint aware correction
        Rate5(k,n) = Rate5(k,n)+log2(det(eye(Sysparam.Ns)+(Rho/Sysparam.Ns)*H(:,:,k)*F3*V5(:,1:Sysparam.Ns,k)*V5(:,1:Sysparam.Ns,k)'*F3'*H(:,:,k)'));
        %% without beam squint compensation
        Rate6(k,n) = Rate6(k,n)+log2(det(eye(Sysparam.Ns)+(Rho/Sysparam.Ns)*H(:,:,k)*F3*Veff5(:,1:Sysparam.Ns,k)*Veff5(:,1:Sysparam.Ns,k)'*F3'*H(:,:,k)'));
    end  
end
clear idealF 
clear F3
end 
%%
Rate2 = real(Rate2)./Nrlz;
Rate3 = real(Rate3)./Nrlz;
Rate5 = real(Rate5)./Nrlz;
Rate6 = real(Rate6)./Nrlz;
%% 
% 
%minRate =[min(Rate2);min(Rate3);min(Rate5);min(Rate6)] 
% AverageRate1 = sum(real(Rate1))./(Nrlz*Sysparam.K)
% AverageRate2 = sum(real(Rate2))./(Nrlz*Sysparam.K)
% AverageRate3 = sum(real(Rate3))./(Nrlz*Sysparam.K)
%%
set(0, 'defaultlinelinewidth', 3); set(0, 'defaultlinemarkersize', 8);
set(0, 'defaultaxesfontsize', 20); set(0, 'defaulttextfontsize', 20); 
for n = 1:length(NantTX)
figure(n)
f2 = cdfplot(Rate2(:,n));
%cftool( get(f2,'XData'), get(f2,'YData') )
set(f2,'LineWidth',4,'LineStyle','-','Marker','none');
hold on
f3 = cdfplot(Rate3(:,n));
%cftool( get(f2,'XData'), get(f2,'YData') )
set(f3,'LineWidth',4,'LineStyle',':','Marker','none');
hold on
% f5 = cdfplot(Rate5(:,n));
% set(f5,'LineWidth',4,'LineStyle','-','Marker','none');
% f6 = cdfplot(Rate6(:,n));
% set(f6,'LineWidth',4,'LineStyle','-','Marker','none');
% hold on
legend('Proposed approach','Tan');
%legend('Proposed approach','BSA','Without compensation');
xlabel('Achievable rate (bit/s/Hz)');
ylabel('CDF');
title("N_t = " + NantTX(n));
end 
%%
minRate =[min(Rate2);min(Rate3);min(Rate5);min(Rate6)] 
maxRate =[max(Rate2);max(Rate3);max(Rate5);max(Rate6)] 
avgRate =[mean(Rate2,1);mean(Rate3,1);mean(Rate5,1);mean(Rate6,1)]