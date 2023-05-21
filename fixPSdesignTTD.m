%input:
% Sysparam and Simparam
%output:
% - PS precoder  F1: dimensions
% - TTD precoder F2k: dimensions ()
% PS values x
%
function [F1,F2k,x,TTD_Tan] = fixPSdesignTTD(Sysparam, Simparam)
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

%% fixed PS precoder
 %% generate phase values
 % MN = Nt
    for l = 1:L
        for n = 1:N
            for m = 1:M
                x(l,n,m) = -(n-1)*sin(phi(l)); 
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

%% generate time-delay values
    for l = 1:L
        if sin(phi(l)) > 0
            for m = 1:M
                TTD_Tan(l,m) = m*(N*sin(phi(l))/2)/fc;
                if TTD_Tan(l,m) >= tmax
                    TTD_Tan(l,m) = tmax;
                end
            end
        else
            for m = 1:M
                TTD_Tan(l,m) = M*abs(N*sin(phi(l)))/(2*fc) + m*(N*sin(phi(1))/2)/fc;
                if TTD_Tan(l,m) >= tmax
                    TTD_Tan(l,m) = tmax;
                end
            end 
        end
    end
    %% TTD precoder
    %% TTD precoder
    F2k = zeros(M*L,L,K);
    for k = 1:K
      for l = 1:L
          F2k((l-1)*M+1:l*M,l,k) =exp(-1j*2*pi*F(k).*TTD_Tan(l,:)); 
      end 
    end 
end 