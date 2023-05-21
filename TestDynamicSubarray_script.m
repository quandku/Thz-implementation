clc;
clear all;
close all;

%% system parameters
Sysparam = struct();

Sysparam.Nt = 256;
Sysparam.Nr = 4;
Sysparam.Nrf=4;
Sysparam.Ns = 4;
Sysparam.fc = 300*1e9;%carrier frequency 
Sysparam.f = 30*1e9; %bandwidth 
Sysparam.K = 129;  %number of sub-carriers
Sysparam.M = 16; %number of TTDs for each RF chain
Sysparam.N = Sysparam.Nt/Sysparam.M; %number of PSs for each TTD
Sysparam.L = 4 ; % number of paths equals to number of RF chains
Sysparam.tmax = 250e-12;% maximum time-delay value
% compute discrete time-delay values
Sysparam.Q = 16; %quantization level 
Sysparam.T = Sysparam.tmax*(0:1:Sysparam.Q-1)/(Sysparam.Q-1);
%% Sub-carrier frequencies
for k =1:Sysparam.K 
    Sysparam.F(k,1) = Sysparam.fc +(Sysparam.f/Sysparam.K)*(k-1-(Sysparam.K-1)/2); % sub-carrier frequency
end 
% relative frequencies
Sysparam.xi = Sysparam.F/Sysparam.fc; 

%% simulation parameters
Simparam = struct();

%Simparam.phi = pi*rand(Sysparam.L,1)-0.5*pi; %AoD;
Simparam.phi = 0.45*pi*ones(Sysparam.L,1)%AoD;
%Simparam.phi = [0.62;0.62;0.62;0.62];
%Simparam.phi = 0.2*pi*ones(Sysparam.L,1); %AoD;
Simparam.psi = pi*rand(Sysparam.L,1)-0.5*pi; %AoA;
Simparam.P = -20:5:15;% transmit power in dB
Simparam.Rho = 10.^(Simparam.P/10); %transmit power
Simparam.Niter = 50; %channel realizations

%generate Fixed TTD matrix  
  F2k = zeros(Sysparam.Q*Sysparam.L,Sysparam.L,Sysparam.K);
    for k = 1:Sysparam.K
        for l = 1:Sysparam.L
            F2k((l-1)*Sysparam.Q+1:l*Sysparam.Q,l,k) = exp(-1j*2*pi*Sysparam.F(k).*Sysparam.T)';
        end 
    end 

%Generate wideband channel
H = WidebandChannel(Sysparam,Simparam);
for k = 1:Sysparam.K
   %[U(:,:,k),D(:,:,k),V(:,:,k)] = svd(H(:,:,k),'econ');
   [~,~,V(:,:,k)] = svd(H(:,:,k),'econ'); %find the optimal digital precoder
end

%Initiate the swicth matrix
S = zeros(Sysparam.Nt, Sysparam.L*Sysparam.Q);
index = randi([1,Sysparam.L*Sysparam.Q],1,Sysparam.Nt);
for n = 1:Sysparam.Nt
    S(n,index(n)) = 1;
end
%initiate digital precoder
% design digital precoder based on the initiated S
        for k = 1:Sysparam.K
            [Uhat(:,:,k),~,Vhat(:,:,k)] = svd(V(:,:,k)'*S*F2k(:,:,k),'econ'); 
            D(:,:,k) = Vhat(:,:,k)*Uhat(:,:,k)';
        end
S1 = S;

 %% Design S and D
tic
 Nloop = 50;
    for iter = 1:Nloop
        % update switch
        for n = 1:Sysparam.Nt
            Row(:,n) = zeros(Sysparam.L*Sysparam.Q,1);
            for k=1:Sysparam.K
                 Row(:,n) = Row(:,n) + (-2*real(F2k(:,:,k)*D(:,:,k)*V(n,:,k)') + diag(F2k(:,:,k)*D(:,:,k)*D(:,:,k)'*F2k(:,:,k)')); 
            end
            [~,imax(n)] = max(Row(:,n));
            S(n,:) =zeros(1,Sysparam.L*Sysparam.Q);
            S(n,imax(n)) = 1;
        end
        % design digital precoder based on the updated S
        for k = 1:Sysparam.K
            [Uhat(:,:,k),~,Vhat(:,:,k)] = svd(V(:,:,k)'*S*F2k(:,:,k),'econ'); 
            D(:,:,k) = Vhat(:,:,k)*Uhat(:,:,k)';
            Loss1(iter,k) = norm(V(:,:,k)-S*F2k(:,:,k)*D(:,:,k),"fro");
        end
        Loss(iter) = sum(Loss1(iter,:));
    end 
toc
    S2 =S;
%%
plot(1:Nloop,Loss)
