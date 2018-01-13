function SINR_output = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, MF_ZF)

    SINR_temp = zeros(K, L);
    
    H_est = zeros(M, K, L);
    for i = 1:L
        for k = 1:K
            H_est(:,k,i) = H(:,k,i,i);                                  % effective channel signal
            for i1 = 1:L
                for k1 = 1:K
                    if P(i1,k1) == P(i,k) && i1~=i && P(i,k)<100 && P(i,k)~=-1
                        H_est(:,k,i) = H_est(:,k,i) + H(:,k1,i,i1);     % inter-cell interference of channel estimation
                    end
                    if  P(i1,k1) == P(i,k)&& P(i,k)>100 && (abs(i1-i)==7||abs(i1-i)==14)
                             H_est(:,k,i) = H_est(:,k,i) + H(:,k1,i,i1);
                    end
                end
            end
            H_est(:,k,i) = H_est(:,k,i) + 1/sqrt(rho_pilot)*random('norm', 0, 1, M, 1);
        end
    end
    
    for l = 1:L
        if MF_ZF == 1
            A = H_est(:,:,l)';
        else
            A = pinv(H_est(:,:,l));
        end

        for k = 1:K
            ICI1 = 0;                                                       % Intra-cell interference
            for k1 = 1:K
                if k1~=k
                    ICI1 = ICI1 + rho_ul*norm(A(k,:)*H(:,k1,l,l))^2;
                end
            end
            ICI2 = 0;                                                       % Inter-cell interference
            for j1 = 1:L
                for k1 = 1:K
                    if j1~=l
                        ICI2 = ICI2 + rho_ul*norm(A(k,:)*H(:,k1,l,j1))^2;
                    end
                end
            end
            ICIN = norm(A(k,:))^2;                                          % equivalent noise
            SINR_temp(k, l) = rho_ul*norm(A(k,:)*H(:,k,l,l))^2/(ICI1 + ICI2 + ICIN);
        end
    end
    SINR_output = SINR_temp;
end
