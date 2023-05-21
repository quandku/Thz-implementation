%generate wideband multiple channel
%phi AoD
%psi AoA
%NT number of transmit antennas
%NR number of received antennas
%F subcarrier frequencies Sysparam,Simparam
function [H] = WidebandChannel(Sysparam,Simparam)
    F = Sysparam.F;
    L = Sysparam.L;
    NT = Sysparam.Nt;
    NR = Sysparam.Nr;
    xi = Sysparam.xi;
    phi = Simparam.phi;
    psi = Simparam.psi;
    K = Sysparam.K; %number of subcarrier frequencies
    H = zeros(NR,NT,K);
    Td = (20e-9)*randn(L,1);%path delay
    %%
    d = 50; % transmit distance (50 m);
    kabs = 0.0033; %(m^-1) absorption coefficient
    % for k =1:K
    % h_gain(k) = 3*1e8*exp(-0.5*kabs*d)/(4*pi*(F(k))*d); 
    % Alpha(:,k) = ((randn(L,1)+randn(L,1)*1j)./(sqrt(2))).*h_gain(k); %channel path gain
    % end
    Alpha = (randn(L,K)+randn(L,K)*1j)./(sqrt(2)); %channel path gain
    %Alpha = randn(L,1)+randn(L,1)*1j; %channel path gain
    for k = 1:K
        for l = 1:L
            H(:,:,k) = H(:,:,k) + Alpha(l,k)*exp(-j*2*pi*Td(l)*F(k))*ArrResponseGenerateH(NR,psi(l),xi(k))*ArrResponseGenerateH(NT,phi(l),xi(k))';
        end 
        H(:,:,k) = H(:,:,k)./sqrt(L); 
    end 
end 
