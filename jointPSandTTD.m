%input:
% System parameters and simulation parameters
%output:
% - PS precoder  F1: dimensions of
% - TTD precoder F2k: dimensions of 
% PS values x
%
function [F1,F2k,x,t] = jointPSandTTD(Sysparam, Simparam)
%Sysparam : Nt,M,N,L,K,F,fc,tmax
%Simparam : phi
%%
Nt = Sysparam.Nt;
M = Sysparam.M;
N = Sysparam.N;
L = Sysparam.L;
K = Sysparam.K;
F = Sysparam.F;
fc = Sysparam.fc;
tmax = Sysparam.tmax;
phi = Simparam.phi; 
VarTheta = 2*fc*tmax;
  %% Generate time-delay values
    for l= 1:L
        for m = 1:M
            if 0 <= sin(phi(l)) && sin(phi(l)) <= (4*fc*tmax)/((2*m-1)*N-1) 
                t(l,m) = ((2*m-1)*N-1)*sin(phi(l))/(4*fc);
            elseif -(4*fc*tmax)/((2*m-1)*N-1) < sin(phi(l))&& sin(phi(l)) < 0
                t(l,m) = ((2*m-1)*N-1)*sin(phi(l))/(4*fc) + tmax;
            elseif sin(phi(l)) > (4*fc*tmax)/((2*m-1)*N-1)
                t(l,m) = tmax;
            elseif sin(phi(l)) < - (4*fc*tmax)/((2*m-1)*N-1)
                t(l,m) = 0;
            end 
        end
        
    end 
    %% TTD precoder
    F2k = zeros(M*L,L,K);
    for k = 1:K
      for l = 1:L
          F2k((l-1)*M+1:l*M,l,k) =exp(-1j*2*pi*F(k).*t(l,:)); 
      end 
    end  
   %% generate phase values
    for l = 1:L
    for n = 1:N
        for m = 1:M
            %when the spatial direction > 0
            if 0 <= sin(phi(l)) && sin(phi(l)) <= (4*fc*tmax)/((2*m-1)*N-1) 
                x(l,n,m) = (N-2*n+1)*sin(phi(l))/2;
            elseif sin(phi(l)) > (4*fc*tmax)/((2*m-1)*N-1)
                x(l,n,m) = VarTheta - ((m-1)*N+n-1)*sin(phi(l));
            %when the spatial direction <0
            elseif -(4*fc*tmax)/((2*m-1)*N-1) < sin(phi(l)) && sin(phi(l)) < 0
                x(l,n,m) = (N-2*n+1)*sin(phi(l))/2;
            elseif  sin(phi(l)) <= -(4*fc*tmax)/((2*m-1)*N-1) 
                x(l,n,m) = VarTheta - ((m-1)*N+n-1)*sin(phi(l));
            end 
        end 
    end 
    end 
    %% PS precoder
    F1 = zeros(Nt,M*L);
    for l = 1:L
        G(:,:,l) = zeros(Nt,M);
        for m = 1:M
            G((m-1)*N+1:m*N,m,l) = exp(1j*pi.*x(l,:,m)); 
        end 
        F1(:,(l-1)*M+1:l*M) = G(:,:,l)/sqrt(Nt); 
    end 
end 