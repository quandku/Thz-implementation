%input: System parameters and simulation parameters
%output: S - Switch network, F2k - TTD precoder, and D - Digital precoder
function[S,F2k,D,Loss] = DynamicSubarray(Sysparam,Simparam,H,Nloop)
    %% generate the fixed TTD precoder
    F2k = zeros(Sysparam.Q*Sysparam.L,Sysparam.L,Sysparam.K);
    for k = 1:Sysparam.K
        for l = 1:Sysparam.L
            F2k((l-1)*Sysparam.Q+1:l*Sysparam.Q,l,k) = exp(-1j*2*pi*Sysparam.F(k).*Sysparam.T)'./sqrt(Sysparam.Nt);
        end 
    end
    %% Find the optimal fully-digital precoder
    for k = 1:Sysparam.K
        %[U(:,:,k),D(:,:,k),V(:,:,k)] = svd(H(:,:,k),'econ');
        [~,~,V(:,:,k)] = svd(H(:,:,k),'econ'); %find the optimal digital precoder
    end
    %% Design S and D

    % Initiate the swicth matrix
    S = zeros(Sysparam.Nt, Sysparam.L*Sysparam.Q);
    index = randi([1,Sysparam.L*Sysparam.Q],1,Sysparam.Nt);
    for n = 1:Sysparam.Nt
        S(n,index(n)) = 1;
    end
    
    %Initiate digital precoder
     % design digital precoder based on the initiated S
        for k = 1:Sysparam.K
            [Uhat(:,:,k),~,Vhat(:,:,k)] = svd(V(:,:,k)'*S*F2k(:,:,k),'econ'); 
            D(:,:,k) = Vhat(:,:,k)*Uhat(:,:,k)';
        end
    %% update S and D iteratively
       % Nloop = 20;
    for iter = 1:Nloop
        % update switch S
        for n = 1:Sysparam.Nt
            Row(:,n) = zeros(Sysparam.L*Sysparam.Q,1);
            for k=1:Sysparam.K
                 Row(:,n) = Row(:,n) + (-2*real(F2k(:,:,k)*D(:,:,k)*V(n,:,k)') + diag(F2k(:,:,k)*D(:,:,k)*D(:,:,k)'*F2k(:,:,k)')); 
            end
            [~,imax(n)] = max(Row(:,n));
            S(n,:) =zeros(1,Sysparam.L*Sysparam.Q);
            S(n,imax(n)) = 1;
        end
        % update digital precoder D(:,:,k)
        for k = 1:Sysparam.K
            [Uhat(:,:,k),~,Vhat(:,:,k)] = svd(V(:,:,k)'*S*F2k(:,:,k),'econ'); 
            D(:,:,k) = Vhat(:,:,k)*Uhat(:,:,k)';
            Loss1(iter,k) = norm(V(:,:,k)-S*F2k(:,:,k)*D(:,:,k),"fro");
        end
        Loss(iter) = sum(Loss1(iter,:));
   end 
end