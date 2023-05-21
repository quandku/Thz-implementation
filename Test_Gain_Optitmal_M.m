%compare beam squint compensation between different method
clc;
clear all;
close all;

%B = 1e9.*[30,40,60];
B = 1e9.*[30];
%B =1e9.*30; 
g0 = 0.9; 

%% system parameters
Sysparam = struct();

Sysparam.Nt = 720;
Mset =divisors(Sysparam.Nt);
Sysparam.Nr = 4;
Sysparam.Ns = 4;
Sysparam.fc = 3e11;%carrier frequency 
BRatio = B./Sysparam.fc;
Sysparam.K = 129;  %number of sub-carriers
Sysparam.L = 1 ; % number of paths equals to number of RF chains
%%
%% simulation parameters
Simparam = struct();
Simparam.psi = pi*rand(Sysparam.L,1)-0.5*pi; %AoA;
Simparam.Niter = 100; %channel realizations
Simparam.phi = 0.2951*pi*ones(Sysparam.L,1);%AoD;
SpatialDirection = sin(Simparam.phi)
for i = 1:length(B)
Omega(i) = 6*(1-g0)/(pi*(Sysparam.K-1)*BRatio(i)*0.8/(4*Sysparam.K))^2;
Mtemp(i) = Sysparam.Nt/sqrt(1+Omega(i));    
mIndex = find(Mset > Mtemp(i));
Mdivisor = divisors(Sysparam.Nt);
MfeasibleSet = Mset(mIndex);  
Mstar(i) = Mset(mIndex(1)); %number of TTDs for each RF chain
Msub(i) = Mset(mIndex(1)-1);
Msup(i) = Mset(mIndex(1)+1);
end
Mtemp; 
Msub;
Mstar;
Msup;
Mcon = [Mstar;Msub];
% compute the poer 
p1 = 100; 
p2 = 20; 


%% 
for i = 1:length(B)
Sysparam.f = B(i); %bandwidth
    for iIter = 1:length(Mdivisor)
    %Sysparam.M = MfeasibleSet(iIter,i)   
    Sysparam.M = Mdivisor(iIter);
    Sysparam.N = Sysparam.Nt/Sysparam.M; %number of PSs for each TTD
    Time(i)= max(sin(Simparam.phi))*((2*Sysparam.M-1)*Sysparam.Nt-Sysparam.M)/(4*Sysparam.M*Sysparam.fc);
    Sysparam.tmax = max(sin(Simparam.phi))*((2*Sysparam.M-1)*Sysparam.Nt-Sysparam.M)/(4*Sysparam.M*Sysparam.fc);
%% Sub-carrier frequencies
for k =1:Sysparam.K 
    Sysparam.F(k,1) = Sysparam.fc +(Sysparam.f/Sysparam.K)*(k-1-(Sysparam.K-1)/2); % sub-carrier frequency
end 
% relative frequencies
Sysparam.xi = Sysparam.F/Sysparam.fc; 
%Generate wideband channel
[PS,TTD,x,t] = jointPSandTTD(Sysparam, Simparam); % joint delay and phase precoding
for l = 1:Sysparam.L
    a(:,l) = arrayResponse(sin(Simparam.phi(l)),Sysparam.Nt);
    for k = 1:Sysparam.K
         b(:,k,l) = arrayResponse(Sysparam.xi(k)*sin(Simparam.phi(l)),Sysparam.Nt);
         arrayGain(k,l,iIter,i) = abs(a(:,l)'*b(:,k,l));
         Compensated1(k,l,iIter,i) = abs((PS*TTD(:,l,k))'*b(:,k,l));
    end
    minGain(l,iIter,i) = min(Compensated1(:,l,iIter,i));
end
clear a b ;
end 
end 
%%
set(0, 'defaultlinelinewidth', 3); set(0, 'defaultlinemarkersize', 8);
set(0, 'defaultaxesfontsize', 15); set(0, 'defaulttextfontsize', 15); 
plot(Mdivisor,minGain)
xline(Msub,'--b','Linewidth',3)
xline(Mstar,'-.r','Linewidth',3)
yline(g0,':k','Linewidth',3)
grid on
xlim([1,Sysparam.Nt])
xlabel('M')
ylabel('Minimum Array Gain')
legend('','M = M^*','M = M_{sub}')
%%

% plot(B./1e9,Mstar,'-*')
% grid on;
% xlabel('Bandwidth B (GHz)')
% ylabel('Minimum number of TTD per RF chain M^*')
% %%
% close all;
% figure(1)
% f1 = cdfplot(Compensated1(:,:,1,1));
% set(f1,'LineWidth',4,'LineStyle',':','Marker','none');
% hold on
% f2 = cdfplot(Compensated1(:,:,2,1));
% set(f2,'LineWidth',4,'LineStyle','--','Marker','none');
% hold on
% f4 = cdfplot(ones(Sysparam.K,1));
% set(f4,'LineWidth',4,'LineStyle','-','Marker','none');
% legend('Proposed approach M = 60','Proposed approach M = 48','Ideal analog precoder')
% xlabel('Array gain')
% ylabel('CDF')
% title('')
% %%
% close all;
% figure(1)
% f1 = cdfplot(Compensated1(:,:,1,1));
% set(f1,'LineWidth',4,'LineStyle',':','Marker','none');
% hold on
% f2 = cdfplot(Compensated1(:,:,2,1));
% set(f2,'LineWidth',4,'LineStyle','--','Marker','none');
% hold on
% f4 = cdfplot(ones(Sysparam.K,1));
% set(f4,'LineWidth',4,'LineStyle','-','Marker','none');
% legend('Proposed approach M = 60','Proposed approach M = 48','Ideal analog precoder')
% xlabel('Array gain')
% ylabel('CDF')
% title('')
