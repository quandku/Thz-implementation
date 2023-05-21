%input:
% - phi: spatial directions
% - fc:
%output: 
function [F1,F2k,x,t] = jointTTDandPS(Sysparam, Simparam)
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
%% 
VarTheta = 2*fc*tmax;
%generate time-delay values
    for l= 1:L
        for m = 1:M
            if 0 <= sin(phi(l)) <= (4*fc*tmax)/((2*m-1)*N-1) 
                t(l,m) = ((2*m-1)*N-1)*sin(phi(l))/(4*fc);
            else
                t(l,m) = tmax;
            end 
        end
    end 
    %generate phase values
    for l = 1:L
    for n = 1:N
        for m = 1:M
            if 0 <= sin(phi(l)) <= (4*fc*tmax)/((2*m-1)*N-1) 
                x(l,n,m) = (N-2*n+1)*sin(phi(l))/2;
            else
                x(l,n,m) = VarTheta - ((m-1)*N+n-1)*sin(phi(l));
            end 
        end 
    end 
    end 
    %Phase-precoder
    F1 = zeros(Nt,M*L);
    for l = 1:L
        G(:,:,l) = zeros(Nt,M);
        for m = 1:M
            G((m-1)*N+1:m*N,m,l) = exp(1j*pi.*x(l,:,m)); 
        end 
        F1(:,(l-1)*M+1:l*M) = G(:,:,l)/sqrt(Nt); 
    end 
    %time-delay precoder
    F2k = zeros(M*L,L,K);
    for k = 1:K
      for l = 1:L
          F2k((l-1)*M+1:l*M,l,k) =exp(-1j*2*pi*F(k).*t(l,:)); 
      end 
    end  

end 