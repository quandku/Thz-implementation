clc;
clear all;
%%
Simparam = struct();

Simparam.M = 256;
Simparam.N = 4;
Simparam.phi1 = pi*rand(Simparam.N,1)-0.5*pi;
Simparam.phi2 = pi*rand(Simparam.N,1)-0.5*pi;
%% construct maxtrice in array manifold
for n =1:Simparam.N
    U(:,n) = exp(j*pi*(0:1:Simparam.M-1)*sin(Simparam.phi1(n)))'/sqrt(Simparam.M);
    V(:,n) = exp(j*pi*(0:1:Simparam.M-1)*sin(Simparam.phi2(n)))'/sqrt(Simparam.M);
end 
%%
A = V'*U*U'*V; 
SigmaA = eigs(A,1)