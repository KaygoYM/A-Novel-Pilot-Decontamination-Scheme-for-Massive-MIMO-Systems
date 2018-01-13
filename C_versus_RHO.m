clear;

clc; 
tic;

%% system parameters
%rho_ul_test = [10^(5/10) 10^(10/10) 10^(15/10) 10^(20/10) 10^(25/10)];                                                        % number of BS antennas M
%rho_ul_test = [10^(2/10) 10^(5/10) 10^(8/10) 10^(11/10) 10^(14/10) 10^(17/10) 10^(20/10)];
rho_ul_test = [10^(2/10) 10^(4/10) 10^(6/10) 10^(8/10) 10^(10/10) 10^(12/10) 10^(14/10) 10^(16/10) 10^(18/10) 10^(20/10)];
M = 512; 
L = 19;
K = 10;
S = 15;
Gamma0 = 0.2;
lambda=0.1;
mu2=0.1;
J=floor(L*K*mu2);
%Test_num = 100;
Test_num = 10;
Pt = 1;
P_n = 1;
R_Cell = 1000;
r_Min = 100;
alpha = 3.8;
sigma_shadow = 8;


 Bw = 1;



UL_SINR_CS_MF = zeros(K, L, length(rho_ul_test), Test_num);
UL_SINR_CSsoft_MF = zeros(K, L, length(rho_ul_test), Test_num);
UL_SINR_WGCPA_MF = zeros(K, L, length(rho_ul_test), Test_num);
UL_SINR_WGCPAIM_MF = zeros(K, L, length(rho_ul_test), Test_num);
UL_SINR_WGCPAsoft_MF = zeros(K, L, length(rho_ul_test), Test_num);




DL_SINR_CS_MF = zeros(K, L, length(rho_ul_test), Test_num);
DL_SINR_CSsoft_MF = zeros(K, L, length(rho_ul_test), Test_num);
DL_SINR_WGCPA_MF = zeros(K, L, length(rho_ul_test), Test_num);
DL_SINR_WGCPAsoft_MF = zeros(K, L, length(rho_ul_test), Test_num);
DL_SINR_GCPA_MF = zeros(K, L, length(rho_ul_test), Test_num);

%% Simulation
for i_M = 1:length(rho_ul_test)
    
    rho_pilot = rho_ul_test(i_M);
    display(rho_pilot);
    rho_dl =rho_pilot;                                                    
    rho_ul =rho_pilot;
    for i_test = 1:Test_num
        display(i_test);
        
        %% Generate Channel vector
        [H Beta] = F_H_Generate(M, L, K, R_Cell, r_Min, sigma_shadow, alpha);

        %% Random Pilot Assignment: k-th pilot for k-th user
        P = zeros(L, K);
        for i = 1:L                                                         % random pilot assignment
            temp = randperm(S);
            P(i,:) = temp(1:K);
        end
        
        UL_SINR_CS_MF(:,:,i_M,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);
        
        
  %      DL_SINR_CS_MF(:,:,i_M,i_test) = F_DLSINR_CAL(L, K, M, H, P, rho_pilot, rho_dl, 2);
        
        
         %% Random Pilot Assignment For Certen_SOFT
      
        [P eu]= F_CS_SOFT(L, K, S, Beta,lambda);
%         for i = 1:L                                                         % random pilot assignment
%             temp = randperm(S);
%             P(i,:) = temp(1:K);
%         end
        
        UL_SINR_CSsoft_MF(:,:,i_M,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);
               
    %    DL_SINR_CSsoft_MF(:,:,i_M,i_test) = F_DLSINR_CAL(L, K, M, H, P, rho_pilot, rho_dl, 2);

        
        
        
        %% Weighted Graph Coloring Based Pilot Assignment WGC-PA
        P = F_WGCPA_Pilot(L, K, S, Beta);
        
        UL_SINR_WGCPA_MF(:,:,i_M,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);

 %       DL_SINR_WGCPA_MF(:,:,i_M,i_test) = F_DLSINR_CAL(L, K, M, H, P, rho_pilot, rho_dl, 2);
     
        
      %% Weighted Graph Coloring Based Pilot Assignment For Certen_SOFT
        [P eu]= F_WGCPA_soft_Pilot(L, K, S, Beta,lambda);
        
        UL_SINR_WGCPAsoft_MF(:,:,i_M,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);

    %    DL_SINR_WGCPAsoft_MF(:,:,i_M,i_test) = F_DLSINR_CAL(L, K, M, H, P, rho_pilot, rho_dl, 2);
        
      %% IM_Weighted Graph Coloring Based Pilot Assignment WGC-PA
        P2 = F_WGCPA_Pilot_IM(L, K, S, Beta, J);
        
        UL_SINR_WGCPAIM_MF(:,:,i_M,i_test) = F_ULSINR_CAL(L, K, M, H, P2, rho_pilot, rho_ul, 1);
        
    end        
end


UL_R_CS_MF = zeros(length(rho_ul_test), 1);
UL_R_CSsoft_MF = zeros(length(rho_ul_test), 1);
UL_R_WGCPA_MF = zeros(length(rho_ul_test), 1);
UL_R_WGCPAIM_MF = zeros(length(rho_ul_test), 1);
UL_R_WGCPAsoft_MF = zeros(length(rho_ul_test), 1);


DL_R_CS_MF = zeros(length(rho_ul_test), 1);
DL_R_CSsoft_MF = zeros(length(rho_ul_test), 1);
DL_R_WGCPA_MF = zeros(length(rho_ul_test), 1);
DL_R_WGCPAsoft_MF = zeros(length(rho_ul_test), 1);
DL_R_GCPA_MF = zeros(length(rho_ul_test), 1);

for i_M = 1:length(rho_ul_test)
    for i_test = 1:Test_num
        for i = 1:L
            for k = 1:K
                UL_R_CS_MF(i_M) = UL_R_CS_MF(i_M) + log2(1+UL_SINR_CS_MF(k,i,i_M,i_test))/Test_num/L/K;
                UL_R_CSsoft_MF(i_M) = UL_R_CSsoft_MF(i_M) + log2(1+UL_SINR_CSsoft_MF(k,i,i_M,i_test))/Test_num/L/K;
                UL_R_WGCPA_MF(i_M) = UL_R_WGCPA_MF(i_M) + log2(1+UL_SINR_WGCPA_MF(k,i,i_M,i_test))/Test_num/L/K;
                UL_R_WGCPAIM_MF(i_M) = UL_R_WGCPAIM_MF(i_M) + log2(1+UL_SINR_WGCPAIM_MF(k,i,i_M,i_test))/Test_num/L/K;
                UL_R_WGCPAsoft_MF(i_M) = UL_R_WGCPAsoft_MF(i_M) + log2(1+UL_SINR_WGCPAsoft_MF(k,i,i_M,i_test))/Test_num/L/K;
               

%                 DL_R_CS_MF(i_M) = DL_R_CS_MF(i_M) + log2(1+DL_SINR_CS_MF(k,i,i_M,i_test))/Test_num/L/K;
%                 DL_R_CSsoft_MF(i_M) = DL_R_CSsoft_MF(i_M) + log2(1+DL_SINR_CSsoft_MF(k,i,i_M,i_test))/Test_num/L/K;
%                 DL_R_WGCPA_MF(i_M) = DL_R_WGCPA_MF(i_M) + log2(1+DL_SINR_WGCPA_MF(k,i,i_M,i_test))/Test_num/L/K;
%                 DL_R_WGCPAsoft_MF(i_M) = DL_R_WGCPAsoft_MF(i_M) + log2(1+DL_SINR_WGCPAsoft_MF(k,i,i_M,i_test))/Test_num/L/K;
%                 DL_R_GCPA_MF(i_M) = DL_R_GCPA_MF(i_M) + log2(1+DL_SINR_GCPA_MF(k,i,i_M,i_test))/Test_num/L/K;
              
            end
        end
    end
end


lw = 1.0;
Tau = 20;
%rho_ul_test=[5 10 15 20 25];
%rho_ul_test=[2 5 8 11 14 17 20];
rho_ul_test=[2 4 6 8 10 12 14 16 18 20];
%% Uplink
figure;

plot(rho_ul_test, UL_R_CS_MF, 'bo-', 'Linewidth', lw);
hold on;
plot(rho_ul_test, UL_R_CSsoft_MF, 'k>-', 'Linewidth', lw);
hold on;
plot(rho_ul_test, UL_R_WGCPA_MF, 'r^-', 'Linewidth', lw);
hold on;
plot(rho_ul_test, UL_R_WGCPAIM_MF, 'mx-', 'Linewidth', lw);
hold on;
plot(rho_ul_test, UL_R_WGCPAsoft_MF, 'cs-', 'Linewidth', lw);
hold on;


grid on;
legend('随机导频分配','软导频分配','加权图染色导频分配','改进加权图染色导频分配','软导频-加权图染色导频分配', 'Location', 'Best');
xlabel('发送功率 (dB)');
ylabel('用户上行平均容量 (bps/Hz)');
hold off;




%% Downlink
% figure(2); 
% 
% plot(rho_ul_test, Bw*T_ul/T_0*DL_R_CS_MF, 'bo-', 'Linewidth', lw);
% hold on;
% plot(rho_ul_test, Bw*T_ul/T_0*DL_R_CSsoft_MF, 'k>-', 'Linewidth', lw);
% hold on;
% plot(rho_ul_test, Bw*T_ul/T_0*DL_R_WGCPA_MF, 'r^-', 'Linewidth', lw);
% hold on;
% plot(rho_ul_test, Bw*T_ul/T_0*DL_R_WGCPAsoft_MF, 'cs-', 'Linewidth', lw);
% hold on;
% plot(rho_ul_test, Bw*T_ul/T_0*DL_R_GCPA_MF, 'gd-', 'Linewidth', lw);
% hold on;
% 
% grid on;
% legend('随机导频分配', '随机导频分配soft','加权图染色导频分配','加权图染色导频分配soft','图染色导频分配',  'Location', 'Best');
% xlabel('发送功率 (dB)');
% ylabel('用户下行可达速率 (Mbits/s/user)');
% hold off;

toc;















