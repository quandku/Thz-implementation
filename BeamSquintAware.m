function [F1,F2k,D] = BeamSquintAware(Sysparam,Simparam,H,Nloop)
%%Construct frequency-independent Analog precoder and Frequency-Dependent
%%Analog precoders
for l = 1:Sysparam.L
    F_Fre_Ind(:,l) = ArrResponseGenerateH(Sysparam.Nt,Simparam.phi(l),1)/sqrt(Sysparam.Nt);
    for k =1:K
        idealF(:,l,k) = ArrResponseGenerateH(Sysparam.Nt,Simparam.phi(l),Sysparam.xi(k))/sqrt(Sysparam.Nt);
    end 
end 
%% 
end