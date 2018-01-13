function P_output = F_WGCPA_Pilot_IM(L, K, S, Beta,J)
    
    %eu=2;
    eta = zeros(L,K,L,K);
    gama=zeros(L,K);
    P = zeros(L, K);
    


    beta_2=zeros(L,K);
    
    for l=1:L
       for k=1:K
        beta_2(l,k)=Beta(k,l,l)^2;
       end
    end
    
    for j=1:J
        a=min(min(beta_2));
        [j0,k0]=find(beta_2==a);
        P(j0,k0)=-1;
        beta_2(j0,k0)=100000;
    end

    for i1 = 1:L                                                            % pilot contamination strength calculation
        for k1 = 1:K
            for i2 = 1:L
                for k2 = 1:K
                    if (i1~=i2)
                        eta(i1,k1,i2,k2) = Beta(k1,i2,i1)^2/Beta(k2,i2,i2)^2+Beta(k2,i1,i2)^2/Beta(k1,i1,i1)^2;
                    end
                end
            end
        end
    end
    
    
%     for i=1:L
%         for k=1:K
%             for ii=1:L
%                 for kk=1:K
%                      gama(i,k)=gama(i,k)+eta(i,k,ii,kk);
%                 end
%             end
%         end
%     end
%     
%     for j=1:J
%         a=max(max(gama));
%         [j0,k0]=find(gama==a);
%         P(j0,k0)=-1;
%         gama(j0,k0)=-100000;
%     end
%     ma_eta=max(reshape(eta,1,L*K*L*K));
%     mi_eta=min(reshape(eta,1,L*K*L*K));
%     rho=mi_eta+

    
 
    eta_max = 0;
    for j_bs = 1:L                                                          % First 2 users pilot assignment
        for k_user = 1:K
            if P(j_bs,k_user)~=-1
             for j1_bs = 1:L
                for k1_user = 1:K
                    if P(j1_bs,k1_user)~=-1
                      if j_bs~=j1_bs && eta(j_bs,k_user,j1_bs,k1_user)>eta_max
                        j_0 = j_bs;
                        k_0 = k_user;
                        j_1 = j1_bs;
                        k_1 = k1_user;
                        eta_max = eta(j_bs,k_user,j1_bs,k1_user);
                      end
                    end
                end
             end
            end
        end
    end
    P(j_0,k_0) = 1;
    P(j_1,k_1) = 2;

    for t = (J+3):L*K   
    %for t = (3+eu*L):L*K  % Pilot assignment for other users
        delta = zeros(L,K);                                                 % UE selection
        delta_max = 0;
        for j = 1:L
            for k = 1:K
                if P(j,k) == 0
                    for j1 = 1:L
                        for k1 = 1:K
                            if j~=j1 && P(j1,k1)>0
                                delta(j,k) = delta(j,k) + eta(j,k,j1,k1);
                            end
                        end
                    end
                    if delta(j,k)>delta_max
                        j_0 = j;
                        k_0 = k;
                        delta_max = delta(j,k);
                    end
                end
            end
        end

        % Pilot selection
        c = zeros(1,S);
        for k = 1:K
            if P(j_0,k)>0
                c(P(j_0,k)) = 1;
            end
        end
        zeta = zeros(1,S);
        for j = 1:L
            for k = 1:K
                if P(j,k)>0
                    zeta(P(j,k)) = zeta(P(j,k)) + eta(j,k,j_0,k_0);
                end
            end
        end
        zeta_min = 1000;
        c0 = 0;
        for k = 1:S
            if zeta(k)<zeta_min && c(k)==0
                zeta_min = zeta(k);
                c0 = k;
            end
        end
        P(j_0,k_0) = c0;
    end
    
    P_output = P;
end