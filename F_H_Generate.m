function [H_output Beta_output] = F_H_Generate(M, L, K, R_Cell, r_Min, sigma_shadow, alpha)

    BaseP = zeros(L,2);                                                         % coords of BS
    BaseP(1,1) = 0;
    BaseP(1,2) = 0;
    for i = 2:7                                                                 % Generate coords of BSs
        BaseP(i,1) = BaseP(1,1) + 2*R_Cell*cos(2*pi*(i-1)/(6));
        BaseP(i,2) = BaseP(1,2) + 2*R_Cell*sin(2*pi*(i-1)/(6));
    end
    
        BaseP(8,1) = BaseP(7,1) + 2*R_Cell*cos(2*pi*(5)/(6));
        BaseP(8,2) = BaseP(7,2) + 2*R_Cell*sin(2*pi*(5)/(6));
        BaseP(9,1) = BaseP(7,1) + 2*R_Cell*cos(2*pi*(6)/(6));
        BaseP(9,2) = BaseP(7,2) + 2*R_Cell*sin(2*pi*(6)/(6));
        BaseP(10,1) = BaseP(7,1) + 2*R_Cell*cos(2*pi*(1)/(6));
        BaseP(10,2) = BaseP(7,2) + 2*R_Cell*sin(2*pi*(1)/(6));
        BaseP(11,1) = BaseP(2,1) + 2*R_Cell*cos(2*pi*(1)/(6));
        BaseP(11,2) = BaseP(2,2) + 2*R_Cell*sin(2*pi*(1)/(6));
        BaseP(12,1) = BaseP(2,1) + 2*R_Cell*cos(2*pi*(2)/(6));
        BaseP(12,2) = BaseP(2,2) + 2*R_Cell*sin(2*pi*(2)/(6));
        BaseP(13,1) = BaseP(3,1) + 2*R_Cell*cos(2*pi*(2)/(6));
        BaseP(13,2) = BaseP(3,2) + 2*R_Cell*sin(2*pi*(2)/(6));
        BaseP(14,1) = BaseP(4,1) + 2*R_Cell*cos(2*pi*(2)/(6));
        BaseP(14,2) = BaseP(4,2) + 2*R_Cell*sin(2*pi*(2)/(6));
        BaseP(15,1) = BaseP(4,1) + 2*R_Cell*cos(2*pi*(3)/(6));
        BaseP(15,2) = BaseP(4,2) + 2*R_Cell*sin(2*pi*(3)/(6));
        BaseP(16,1) = BaseP(4,1) + 2*R_Cell*cos(2*pi*(4)/(6));
        BaseP(16,2) = BaseP(4,2) + 2*R_Cell*sin(2*pi*(4)/(6));
        BaseP(17,1) = BaseP(5,1) + 2*R_Cell*cos(2*pi*(4)/(6));
        BaseP(17,2) = BaseP(5,2) + 2*R_Cell*sin(2*pi*(4)/(6));
        BaseP(18,1) = BaseP(6,1) + 2*R_Cell*cos(2*pi*(4)/(6));
        BaseP(18,2) = BaseP(6,2) + 2*R_Cell*sin(2*pi*(4)/(6));
        BaseP(19,1) = BaseP(6,1) + 2*R_Cell*cos(2*pi*(5)/(6));
        BaseP(19,2) = BaseP(6,2) + 2*R_Cell*sin(2*pi*(5)/(6));
        
        
    
    H = zeros(M,K,L,L);                                                 % channel matrix
    Beta = zeros(K,L,L);                                                % large scale fading coefficients
    for i = 1:L                                                         % k-th user in the j-th cell to i-th BS
        for j = 1:L
            for k = 1:K
                x1 = BaseP(i,1);
                y1 = BaseP(i,2);
%                 x2 = random('unif', 0, 2*R_Cell, 1, 1);
%                 y2 = random('unif', 0, 2*R_Cell, 1, 1);
%                 while ((x2-R_Cell)^2+(y2-R_Cell)^2>R_Cell^2) || ((x2-R_Cell)^2+(y2-R_Cell)^2<r_Min^2)
%                     x2 = random('unif', 0, 2*R_Cell, 1, 1);
%                     y2 = random('unif', 0, 2*R_Cell, 1, 1);
%                 end
                x2 = random('unif', -R_Cell, R_Cell, 1, 1);
                y2 = random('unif', -R_Cell, R_Cell, 1, 1);
                while ((x2)^2+(y2)^2>R_Cell^2) || ((x2)^2+(y2)^2<r_Min^2)
                x2 = random('unif', -R_Cell, R_Cell, 1, 1);
                y2 = random('unif', -R_Cell, R_Cell, 1, 1);
                end
                x2 = x2 + BaseP(j,1);
                y2 = y2 + BaseP(j,2);
                r = sqrt((x1-x2)^2+(y1-y2)^2);                          % distance between UE and BS
                z = 10^(random('norm', 0, 10^(sigma_shadow/20), 1, 1)/10); 
                Beta(k,i,j) = z/((r/R_Cell)^alpha+1);
                H(:,k,i,j) = sqrt(Beta(k,i,j))*1/sqrt(2)*(random('norm', 0, 1, 1, M)+random('norm', 0, 1, 1, M)*1j);
            end
        end
    end
    H_output = H;
    Beta_output = Beta;
end