%input:  - phi =sin(AoD): spatial directions
% - N: number of Tx/Rx antennas
function [ar] = arrayResponse(phi,N)
    ar = exp(j*pi*(0:1:N-1)*phi)'/sqrt(N);
end 