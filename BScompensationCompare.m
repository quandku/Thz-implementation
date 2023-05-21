%compare beam squint compensation between different methods
clc;
clear all;
close all;
%use different values of tmax
TMax = (1e-12)*(200:10:400); 
%% system parameters
Sysparam = struct();

Sysparam.Nt = 256;
Sysparam.Nr = 4;
Sysparam.Ns = 4;
Sysparam.fc = 3000*1e9;%carrier frequency 
Sysparam.f = 300*1e9; %bandwidth 
Sysparam.K = 129;  %number of sub-carriers
Sysparam.M = 16; %number of TTDs for each RF chain
Sysparam.N = Sysparam.Nt/Sysparam.M; %number of PSs for each TTD
Sysparam.L = 1 ; % number of paths equals to number of RF chains
%Sysparam.Nrf=Sysparam.L;
% compute discrete time-delay values
Sysparam.Q = 16; % number of fixed TTD element 
%% Sub-carrier frequencies
for k =1:Sysparam.K 
    Sysparam.F(k,1) = Sysparam.fc +(Sysparam.f/Sysparam.K)*(k-1-(Sysparam.K-1)/2); % sub-carrier frequency
end 
% relative frequencies
Sysparam.xi = Sysparam.F/Sysparam.fc; 

%% simulation parameters
Simparam = struct();
Nloop = 20;
%Simparam.phi = [0.62;0.62;0.62;0.62];
%Simparam.phi = 0.2*pi*ones(Sysparam.L,1); %AoD;
Simparam.psi = pi*rand(Sysparam.L,1)-0.5*pi; %AoA;
Simparam.P = -20:5:15;% transmit power in dB
Simparam.Rho = 10.^(Simparam.P/10); %transmit power
Simparam.Niter = 50; %channel realizations

%Simparam.phi = pi*rand(Sysparam.L,1)-0.5*pi %AoD;
Simparam.phi = 0.2951*pi*ones(Sysparam.L,1)%AoD;
%Simparam.phi = pi*((0.5-1/3)*rand(Sysparam.L,1)+1/3) %AoD;
SpatialDirection = sin(Simparam.phi)

%Generate wideband channel
H = WidebandChannel(Sysparam,Simparam);
%% 
 for i = 1:length(TMax)
     
Sysparam.tmax = TMax(i);% maximum time-delay value
Sysparam.T = Sysparam.tmax*(0:1:Sysparam.Q-1)/(Sysparam.Q-1);
[G1,G2k,x2,t2] = fixPSdesignTTD(Sysparam,Simparam); %Tan design
[PS,TTD,x,t] = jointPSandTTD(Sysparam, Simparam); % joint delay and phase precoding
[S,F2k,D,~] = DynamicSubarray(Sysparam,Simparam,H,Nloop); %Dynamic subarray method
%% change Q
Sysparam.Q = 64; % number of fixed TTD element 
Sysparam.T = Sysparam.tmax*(0:1:Sysparam.Q-1)/(Sysparam.Q-1);
[S1,F2k1,D1,~] = DynamicSubarray(Sysparam,Simparam,H,Nloop); %Dynamic subarray method

Sysparam.Q = 128; % number of fixed TTD element 
Sysparam.T = Sysparam.tmax*(0:1:Sysparam.Q-1)/(Sysparam.Q-1);
[S2,F2k2,D2,~] = DynamicSubarray(Sysparam,Simparam,H,Nloop); %Dynamic subarray method

for l = 1:Sysparam.L
    a(:,l) = arrayResponse(sin(Simparam.phi(l)),Sysparam.Nt);
    for k = 1:Sysparam.K
         b(:,k,l) = arrayResponse(Sysparam.xi(k)*sin(Simparam.phi(l)),Sysparam.Nt);
         arrayGain(k,l) = abs(a(:,l)'*b(:,k,l));
         Compensated1(k,l) = abs((PS*TTD(:,l,k))'*b(:,k,l));
         Compensated2(k,l) = abs((G1*G2k(:,l,k))'*b(:,k,l));          
         Compensated3(k,l) = abs((S*F2k(:,l,k))'*b(:,k,l))/sqrt(Sysparam.Nt);
         Compensated4(k,l) = abs((S1*F2k1(:,l,k))'*b(:,k,l))/sqrt(Sysparam.Nt);
         Compensated5(k,l) = abs((S2*F2k2(:,l,k))'*b(:,k,l))/sqrt(Sysparam.Nt);
    end
end
% plot(1:Sysparam.K,Compensated,'-*',1:Sysparam.K,Compensated2,1:Sysparam.K,arrayGain)
% legend('Proposed method','Tan&Dai design','Conventional hybrid precoding')
%plot(1:Sysparam.K,Compensated1,'-*',1:Sysparam.K,Compensated2,1:Sysparam.K,Compensated3,1:Sysparam.K,arrayGain)
%legend('Proposed method','Tan&Dai design','Dynamic sub-array method','Conventional hybrid precoding')

%compute the average array gain
AverageGain1(i) = sum(Compensated1)/Sysparam.K;
AverageGain2(i) = sum(Compensated2)/Sysparam.K;
AverageGain3(i) = sum(Compensated3)/Sysparam.K;
AverageGain4(i) = sum(Compensated4)/Sysparam.K;
AverageGain5(i) = sum(Compensated5)/Sysparam.K;
%% 
min1(i) = min(Compensated1);
Max1(i) = max(Compensated1);
R1(i) = max(Compensated1)/min(Compensated1); 
MaxA1(i) = max(Compensated1)/AverageGain1(i);
minA1(i) = min(Compensated1)/AverageGain1(i);
min2(i) = min(Compensated2);
Max2(i) = max(Compensated2);
R2(i) = max(Compensated2)/min(Compensated2); 
MaxA2(i) = max(Compensated2)/AverageGain2(i);
minA2(i) = min(Compensated2)/AverageGain2(i);
min3(i) = min(Compensated3);
Max3(i) = max(Compensated3);
R3(i) = max(Compensated3)/min(Compensated3); 
MaxA3(i) = max(Compensated3)/AverageGain3(i);
minA3(i) = min(Compensated3)/AverageGain3(i);
 end 
%% 
set(0, 'defaultlinelinewidth', 1); set(0, 'defaultlinemarkersize', 8);
set(0, 'defaultaxesfontsize', 15); set(0, 'defaulttextfontsize', 15); 
plot(TMax*1e+12,AverageGain1,'-*b',TMax*1e+12,AverageGain2,'-ok',TMax*1e+12,sqrt(Sysparam.Nt)*AverageGain3,'-.r',TMax*1e+12,sqrt(Sysparam.Nt)*AverageGain4,'--k',TMax*1e+12,AverageGain5)
grid on;,
legend('Proposed method','fixed PS method','Dynamic subarray method with Q = 32','Dynamic subarray method with Q = 64','Dynamic subarray method with Q = 128')
xlabel('Maximum time delay values t_{max} (ps)')
ylabel('Average array gain')
