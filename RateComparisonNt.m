close all;  
clear variables;
%use different values of tmax
NantN = [32,64,128,256,512];
%NantN = 32;
%% system parameters
Sysparam = struct();

Sysparam.Nr = 4;
Sysparam.Nrf= 4;
Sysparam.Ns = 4;
Sysparam.fc = 300*1e9;%carrier frequency 
Sysparam.f = 30*1e9; %bandwidth 
Sysparam.K = 129;  %number of sub-carriers
Sysparam.M = 32; %number of TTDs for each RF chain

Sysparam.L = 4 ; % number of paths equals to number of RF chains
Sysparam.Q = 32; % number of fixed TTD element
Sysparam.tmax = 300*1e-12;% maximum time-delay value 
Sysparam.T = Sysparam.tmax*(0:1:Sysparam.Q-1)/(Sysparam.Q-1);


%% simulation parameters
Simparam = struct();

%Simparam.phi = pi*rand(Sysparam.L,1)-0.5*pi; %AoD;
Simparam.phi = pi*((0.5-0.2951)*rand(Sysparam.L,1)+0.2951) %AoD;
Simparam.psi = pi*rand(Sysparam.L,1)-0.5*pi; %AoA;

%Simparam.P = -20:5:15;% transmit power in dB
Simparam.P = 5;% transmit power in dB
Simparam.Rho = 10.^(Simparam.P/10); %transmit power
Simparam.Niter = 50; %channel realizations
%% Sub-carrier frequencies
for k =1:Sysparam.K 
    Sysparam.F(k,1) = Sysparam.fc +(Sysparam.f/Sysparam.K)*(k-1-(Sysparam.K-1)/2); % sub-carrier frequency
end 
% relative frequencies
Sysparam.xi = Sysparam.F/Sysparam.fc; 
%% 
%  H = WidebandChannel(Sysparam,Simparam);   
% [S,F2k,D,~] = DynamicSubarray(Sysparam,Simparam,H); %Dynamic subarray method
%% Achievable rate comaprison

Rate1 = zeros(1,length(NantN)); %joint delay and phase precoding
Rate2 = zeros(1,length(NantN)); %fixed PS and design TTD 
Rate3 = zeros(1,length(NantN)); % fixed TTD and design switch network
Rate4 = zeros(1,length(NantN));
Rate5 = zeros(1,length(NantN));
% for k =1:Sysparam.K
%     Gain1(k) = 
% end 
%% 
Nrlz = 50; %channel realizationstic
for n = 1:length(NantN)
    Sysparam.Nt = NantN(n);
    Sysparam.N = Sysparam.Nt/Sysparam.M; %number of PSs for each TTD
    %% without beam squint compensation
    F3 = zeros(Sysparam.Nt,Sysparam.L);
    for k = 1:Sysparam.L
        F3(:,k)= ArrResponseGenerateH(Sysparam.Nt,Simparam.phi(k),1)/sqrt(Sysparam.Nt);
    end
     %% compensate beam squint
    [G1,G2k,~,~] = fixPSdesignTTD(Sysparam,Simparam); %Tan design
    [PS,TTD,~,~] = jointPSandTTD(Sysparam, Simparam); % joint delay and phase precoding
 for l =1:Nrlz
     H = WidebandChannel(Sysparam,Simparam);
%     
%     [S,F2k,D,~] = DynamicSubarray(Sysparam,Simparam,H); %Dynamic subarray method    
     for k = 1:Sysparam.K
        %% Fully digital optimal rate -Rate 1
          [~,~,V(:,:,k)] = svd(H(:,:,k));
          Rate1(n) = Rate1(n)+log2(det(eye(Sysparam.Ns)+(Simparam.Rho/Sysparam.Ns)*H(:,:,k)*V(:,1:Sysparam.Ns,k)*V(:,1:Sysparam.Ns,k)'*H(:,:,k)'));
%         % Joint PS and TTD rate - Rate 2   
%           Heff2(:,:,k) = H(:,:,k)*PS*TTD(:,:,k);
%           [~,~,Veff2(:,:,k)] = svd(Heff2(:,:,k));
%           Rate2(n) = Rate2(n)+log2(det(eye(Sysparam.Ns)+(Simparam.Rho/Sysparam.Ns)*Heff2(:,:,k)*Veff2(:,1:Sysparam.Ns,k)*Veff2(:,1:Sysparam.Ns,k)'*Heff2(:,:,k)'));
%          % Fixed PS and design TTD rate - Rate 3
%          Heff3(:,:,k) = H(:,:,k)*G1*G2k(:,:,k);
%          [~,~,Veff3(:,:,k)] = svd(Heff3(:,:,k));
%          Rate3(n) = Rate3(n)+log2(det(eye(Sysparam.Ns)+(Simparam.Rho/Sysparam.Ns)*Heff3(:,:,k)*Veff3(:,1:Sysparam.Ns,k)*Veff3(:,1:Sysparam.Ns,k)'*Heff3(:,:,k)'));
% %          % Dynamic subarray rate - Rate 4
% %          Heff4(:,:,k) = H(:,:,k)*S*F2k(:,:,k);
% %          [~,~,Veff4(:,:,k)] = svd(Heff4(:,:,k));
% %          Rate4(n) = Rate4(n)+log2(det(eye(Sysparam.Ns)+(Simparam.Rho/Sysparam.Ns)*Heff4(:,:,k)*Veff4(:,1:Sysparam.Ns,k)*Veff4(:,1:Sysparam.Ns,k)'*Heff4(:,:,k)'));    
%          % Without beam squint compensation
%          Heff5(:,:,k) = H(:,:,k)*F3;
%          [~,~,Veff5(:,:,k)] = svd(Heff5(:,:,k));
%          Rate5(n) = Rate5(n)+log2(det(eye(Sysparam.Ns)+(Simparam.Rho/Sysparam.Ns)*Heff5(:,:,k)*Veff5(:,1:Sysparam.Ns,k)*Veff5(:,1:Sysparam.Ns,k)'*Heff5(:,:,k)'));
     end
% Rate1(n) = real(Rate1(n))/(Nrlz*Sysparam.K);
% Rate2(n) = real(Rate2(n))/(Nrlz*Sysparam.K);
% Rate3(n) = real(Rate3(n))/(Nrlz*Sysparam.K);
% Rate4(n) = real(Rate4(n))/(Nrlz*Sysparam.K);
% Rate5(n) = real(Rate5(n))/(Nrlz*Sysparam.K);    
end
end
%%
set(0, 'defaultlinelinewidth', 1); set(0, 'defaultlinemarkersize', 8);
set(0, 'defaultaxesfontsize', 15); set(0, 'defaulttextfontsize', 15); 
plot(TMax,Rate1,'->r',TMax,Rate2,'-*b',TMax,Rate3,'-.k',TMax,Rate4,'-o',TMax,Rate5,'-<')
grid on 
xlabel('SNR (dB)')
ylabel('Achievable rate per subcarrier (bit/s/Hz)')
legend('Optimal fully-digital precoding','Proposed method','fixed PS and design TTD','Dynamic subarray Q = 32','Without beam squint compensation')