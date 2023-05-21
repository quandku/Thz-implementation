function [AR] = ArrResponseGenerateH(N,theta,xi_k)
%theta is AoD/ AoA
%N number of antennas
    n = [1:1:N]';
    AR= exp(-1j*pi*xi_k*sin(theta).*(n-1));
end

