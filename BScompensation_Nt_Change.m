%compare beam squint compensation between different method
clc;
clear all;
close all;
%use different values of tmax
%TMax = (1e-12)*340; %
%NantTX = 2.^(5:1:9);
%NantTX = 32*[1,2,4,5,8,10,12,14,16];
NantTX = 256;

%% system parameters
Sysparam = struct();

%Sysparam.Nt = 256;
Sysparam.Nr = 4;
Sysparam.Ns = 4;
Sysparam.fc = 3e11;%carrier frequency 
Sysparam.f = 30e9; %bandwidth 
Sysparam.K = 129;  %number of sub-carriers
Sysparam.M = 8; %number of TTDs for each RF chain
Sysparam.L = 1 ; % number of paths equals to number of RF chains
% compute discrete time-delay values
Sysparam.Q = 8; % number of fixed TTD element 
Sysparam.tmax = (1e-12)*340;
Sysparam.T = Sysparam.tmax*(0:1:Sysparam.Q-1)/(Sysparam.Q-1);
%% Sub-carrier frequencies
for k =1:Sysparam.K 
    Sysparam.F(k,1) = Sysparam.fc +(Sysparam.f/Sysparam.K)*(k-1-(Sysparam.K-1)/2); % sub-carrier frequency
end 
% relative frequencies
Sysparam.xi = Sysparam.F/Sysparam.fc; 

%% simulation parameters
Simparam = struct();

Simparam.psi = pi*rand(Sysparam.L,1)-0.5*pi; %AoA;
Simparam.P = -20:5:15;% transmit power in dB
Simparam.Rho = 10.^(Simparam.P/10); %transmit power
Simparam.Niter = 50; %channel realizations
Simparam.phi = 0.2951*pi*ones(Sysparam.L,1)%AoD;
SpatialDirection = sin(Simparam.phi)
Nloop = 50;

%% 
 for i = 1:length(NantTX)
     
Sysparam.Nt = NantTX(i);% maximum time-delay value
Sysparam.N = Sysparam.Nt/Sysparam.M; %number of PSs for each TTD
%Generate wideband channel
H = WidebandChannel(Sysparam,Simparam);
[G1,G2k,x2,t2] = fixPSdesignTTD(Sysparam,Simparam); %Tan design
[PS,TTD,x,t] = jointPSandTTD(Sysparam, Simparam); % joint delay and phase precoding
[S,F2k,D,~] = DynamicSubarray(Sysparam,Simparam,H,Nloop); %Dynamic subarray method

for l = 1:Sysparam.L
    a(:,l) = arrayResponse(sin(Simparam.phi(l)),Sysparam.Nt);
    for k = 1:Sysparam.K
         b(:,k,l) = arrayResponse(Sysparam.xi(k)*sin(Simparam.phi(l)),Sysparam.Nt);
         arrayGain(k,l) = abs(a(:,l)'*b(:,k,l));
         Compensated1(k,l) = abs((PS*TTD(:,l,k))'*b(:,k,l));
         Compensated2(k,l) = abs((G1*G2k(:,l,k))'*b(:,k,l));          
         Compensated3(k,l) = abs((S*F2k(:,l,k))'*b(:,k,l));
    end
end
clear a b ;
% plot(1:Sysparam.K,Compensated,'-*',1:Sysparam.K,Compensated2,1:Sysparam.K,arrayGain)
% legend('Proposed method','Tan&Dai design','Conventional hybrid precoding')
%plot(1:Sysparam.K,Compensated1,'-*',1:Sysparam.K,Compensated2,1:Sysparam.K,Compensated3,1:Sysparam.K,arrayGain)
%legend('Proposed method','Tan&Dai design','Dynamic sub-array method','Conventional hybrid precoding')

%compute the average array gain
AverageGain1(i) = sum(Compensated1)/Sysparam.K;
AverageGain2(i) = sum(Compensated2)/Sysparam.K;
AverageGain3(i) = sum(Compensated3)/Sysparam.K;
%% find ratio max-min
min1(i) = min(Compensated1);
Max1(i) = max(Compensated1);
R1(i) = max(Compensated1)/min(Compensated1); 
MaxA1(i) = max(Compensated1)/AverageGain1(i);
minA1(i) = min(Compensated1)/AverageGain1(i);
mMA1(i) = (max(Compensated1)-min(Compensated1))/AverageGain1(i);
min2(i) = min(Compensated2);
Max2(i) = max(Compensated2);
R2(i) = max(Compensated2)/min(Compensated2); 
MaxA2(i) = max(Compensated2)/AverageGain2(i);
minA2(i) = min(Compensated2)/AverageGain2(i);
mMA2(i) = (max(Compensated2)-min(Compensated2))/AverageGain2(i);
min3(i) = min(Compensated3);
Max3(i) = max(Compensated3);
R3(i) = max(Compensated3)/min(Compensated3); 
MaxA3(i) = max(Compensated3)/AverageGain3(i);
minA3(i) = min(Compensated3)/AverageGain3(i);
mMA3(i) = (max(Compensated3)-min(Compensated3))/AverageGain3(i);

 end 
%%
set(0, 'defaultlinelinewidth', 1); set(0, 'defaultlinemarkersize', 8);
set(0, 'defaultaxesfontsize', 15); set(0, 'defaulttextfontsize', 15); 
% figure(1)
% plot(NantTX,AverageGain1,'-*b',NantTX,AverageGain2,'-ok',NantTX,AverageGain3,'->r')
% grid on;
% xlim([0 max(NantTX)])
% legend('Proposed method','fixed PS method','Dynamic subarray method with Q = 32')
% xlabel('Number of transmit antennas N_t')
% ylabel('Average array gain')
% figure(2)
% plot(NantTX,R1,'-*b',NantTX,R2,'-ok',NantTX,R3,'->r')
% grid on;
% xlim([0 max(NantTX)])
% legend('Proposed method','fixed PS method','Dynamic subarray method with Q = 32')
% xlabel('Number of transmit antennas N_t')
% ylabel('Max/Min Ratio')
% figure(3)
% plot(NantTX,MaxA1,'-*b',NantTX,MaxA2,'-ok',NantTX,MaxA3,'->r')
% grid on;
% xlim([0 max(NantTX)])
% legend('Proposed method','fixed PS method','Dynamic subarray method with Q = 32')
% xlabel('Number of transmit antennas N_t')
% ylabel('Max/Average Ratio')
% figure(4)
% plot(NantTX,minA1,'-*b',NantTX,minA2,'-ok',NantTX,minA3,'->r')
% grid on;
% xlim([0 max(NantTX)])
% legend('Proposed method','fixed PS method','Dynamic subarray method with Q = 32')
% xlabel('Number of transmit antennas N_t')
% ylabel('Min/Average Ratio')
figure(5)
plot(1:k,Compensated1,1:k,Compensated2,1:k,Compensated3)
grid on;
xlim([1 129])
legend('Proposed method','fixed PS method','Dynamic subarray method with Q = 32')
xlabel('Subcarrier indices')
ylabel('Array Gain')
% figure(6)
% plot(NantTX,mMA1,'-*b',NantTX,mMA2,'-ok',NantTX,mMA3,'->r')
% grid on;
% xlim([0 max(NantTX)])
% legend('Proposed method','fixed PS method','Dynamic subarray method with Q = 32')
% xlabel('Number of transmit antennas N_t')
% ylabel('Max-Min/Average Ratio')
%%
