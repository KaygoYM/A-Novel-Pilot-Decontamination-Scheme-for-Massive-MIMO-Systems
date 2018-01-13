function SINR_output = F_DLSINR_CAL(L, K, M, H, P, rho_pilot, rho_dl, MF_ZF)

    SINR_temp = zeros(K, L);
    
    H_est = zeros(M, K, L);
    for i = 1:L
        for k = 1:K
            H_est(:,k,i) = H(:,k,i,i);                                  % effective channel signal
            for i1 = 1:L
                for k1 = 1:K
                    if P(i1,k1) == P(i,k) && i1~=i && P(i1,k1)~=-1 && P(i,k)~=-1
                        H_est(:,k,i) = H_est(:,k,i) + H(:,k1,i,i1);     % inter-cell interference of channel estimation
                    end
                end
            end
            H_est(:,k,i) = H_est(:,k,i) + 1/sqrt(rho_pilot)*random('norm', 0, 1, M, 1);
        end
    end
    
    W = zeros(M,K,L);
    for i = 1:L
        if MF_ZF == 1
            W(:,:,i) = conj(H_est(:,:,i));
            W(:,:,i) = W(:,:,i)/sqrt(norm(W(:,:,i),'fro')^2);
        else
            W(:,:,i) = pinv(H_est(:,:,i).');
            W(:,:,i) = W(:,:,i)/sqrt(norm(W(:,:,i),'fro')^2);
        end
    end
    
    H_eff = zeros(K,K,L,L);
    for j = 1:L
        for i = 1:L
            H_eff(:,:,j,i) = H(:,:,j,i).'*W(:,:,j);
        end
    end
    
    for i = 1:L
        for k = 1:K
            ICI1 = 0;                                                       % Intra-cell interference
            for k1 = 1:K
                if k1~=k
                    ICI1 = ICI1 + norm(H_eff(k,k1,i,i))^2;
                end
            end
            ICI2 = 0;                                                       % Inter-cell interference
            for j1 = 1:L
                for k1 = 1:K
                    if j1~=i
                        ICI2 = ICI2 + norm(H_eff(k,k1,j1,i1))^2;
                    end
                end
            end
            ICIN = 1/rho_dl;                                          % equivalent noise
            SINR_temp(k, i) = norm(H_eff(k,k,i,i))^2/(ICI1 + ICI2 + ICIN);
        end
    end
    SINR_output = SINR_temp;
end
