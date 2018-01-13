clear;

clc; 
tic;

%% system parameters
M_test = [32 64 128 256 512 1024 2048];                                                        % number of BS antennas M
L = 19;
K = 10;
S = 15;
Gamma0 = 0.2;
lambda=0.1;
% mu1=0.3;
% eu=floor(K*mu1);
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
rho_ul = 10^(15/10);                                                        % 15 dB
rho_dl = 10^(15/10);                                                        % 15 dB
rho_pilot = 10^(15/10);

 Bw = 1;


UL_SINR_CS_MF = zeros(K, L, length(M_test), Test_num);
UL_SINR_CSsoft_MF = zeros(K, L, length(M_test), Test_num);
UL_SINR_WGCPA_MF = zeros(K, L, length(M_test), Test_num);
UL_SINR_WGCPAIM_MF = zeros(K, L, length(M_test), Test_num);
UL_SINR_WGCPAsoft_MF = zeros(K, L, length(M_test), Test_num);
UL_SINR_GCPA_MF = zeros(K, L, length(M_test), Test_num);



DL_SINR_CS_MF = zeros(K, L, length(M_test), Test_num);
DL_SINR_CSsoft_MF = zeros(K, L, length(M_test), Test_num);
DL_SINR_WGCPA_MF = zeros(K, L, length(M_test), Test_num);
DL_SINR_WGCPAsoft_MF = zeros(K, L, length(M_test), Test_num);
DL_SINR_GCPA_MF = zeros(K, L, length(M_test), Test_num);

%% Simulation
for i_M = 1:length(M_test)
    M = M_test(i_M);
    display(M);
    
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
        
        
   %     DL_SINR_CS_MF(:,:,i_M,i_test) = F_DLSINR_CAL(L, K, M, H, P, rho_pilot, rho_dl, 2);
        
        
         %% Random Pilot Assignment For Certen_SOFT
      
        [P eu1] = F_CS_SOFT(L, K, S, Beta,lambda);
%         for i = 1:L                                                         % random pilot assignment
%             temp = randperm(S);
%             P(i,:) = temp(1:K);
%         end
        
        UL_SINR_CSsoft_MF(:,:,i_M,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);
               
   %     DL_SINR_CSsoft_MF(:,:,i_M,i_test) = F_DLSINR_CAL(L, K, M, H, P, rho_pilot, rho_dl, 2);

        
        
        
        %% Weighted Graph Coloring Based Pilot Assignment WGC-PA
        P = F_WGCPA_Pilot(L, K, S, Beta);
        
        UL_SINR_WGCPA_MF(:,:,i_M,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);

  %      DL_SINR_WGCPA_MF(:,:,i_M,i_test) = F_DLSINR_CAL(L, K, M, H, P, rho_pilot, rho_dl, 2);
     
        
      %% Weighted Graph Coloring Based Pilot Assignment For Certen_SOFT
        [P eu]= F_WGCPA_soft_Pilot(L, K, S, Beta,lambda);
        
        UL_SINR_WGCPAsoft_MF(:,:,i_M,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);

  %      DL_SINR_WGCPAsoft_MF(:,:,i_M,i_test) = F_DLSINR_CAL(L, K, M, H, P, rho_pilot, rho_dl, 2);
   
      %% IM_Weighted Graph Coloring Based Pilot Assignment WGC-PA
        P2 = F_WGCPA_Pilot_IM(L, K, S, Beta, J);
        
        UL_SINR_WGCPAIM_MF(:,:,i_M,i_test) = F_ULSINR_CAL(L, K, M, H, P2, rho_pilot, rho_ul, 1); 
  
    end        
end


UL_R_CS_MF = zeros(length(M_test), 1);
UL_R_CSsoft_MF = zeros(length(M_test), 1);
UL_R_WGCPA_MF = zeros(length(M_test), 1);
UL_R_WGCPAIM_MF = zeros(length(M_test), 1);
UL_R_WGCPAsoft_MF = zeros(length(M_test), 1);

DL_R_CS_MF = zeros(length(M_test), 1);
DL_R_CSsoft_MF = zeros(length(M_test), 1);
DL_R_WGCPA_MF = zeros(length(M_test), 1);
DL_R_WGCPAsoft_MF = zeros(length(M_test), 1);
DL_R_GCPA_MF = zeros(length(M_test), 1);


mu0=0.05;
T_1=1-S/K*mu0;
S_soft=S+eu*7;
T_soft=1-S_soft/K*mu0;
S_WIM=S+J;
T_WIM=1-S_WIM/K*mu0;


for i_M = 1:length(M_test)
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

%% Uplink
figure;

semilogx(M_test, Bw*T_1*UL_R_CS_MF, 'bo-', 'Linewidth', lw);
hold on;
semilogx(M_test, Bw*T_soft*UL_R_CSsoft_MF, 'k>-', 'Linewidth', lw);
hold on;
semilogx(M_test, Bw*T_1*UL_R_WGCPA_MF, 'r^-', 'Linewidth', lw);
hold on;
%semilogx(M_test, Bw*T_WIM*UL_R_WGCPAIM_MF, 'mx-', 'Linewidth', lw);
%hold on;
semilogx(M_test, Bw*T_soft*UL_R_WGCPAsoft_MF, 'cs-', 'Linewidth', lw);
hold on;


grid on;
legend('Random','SPRS', 'WGC-PD','SPRS+WGC-PD', 'Location', 'Best');
xlabel('Number of antennas in BS M');
ylabel('Average uplink achievable rate per user C (bps/Hz)');
hold off;




%% Downlink
% figure(3); 
% 
% semilogx(M_test, Bw*T_ul/T_0*DL_R_CS_MF, 'bo-', 'Linewidth', lw);
% hold on;
% semilogx(M_test, Bw*T_ul/T_0*DL_R_CSsoft_MF, 'k>-', 'Linewidth', lw);
% hold on;
% semilogx(M_test, Bw*T_ul/T_0*DL_R_WGCPA_MF, 'r^-', 'Linewidth', lw);
% hold on;
% semilogx(M_test, Bw*T_ul/T_0*DL_R_WGCPAsoft_MF, 'cs-', 'Linewidth', lw);
% hold on;
% semilogx(M_test, Bw*T_ul/T_0*DL_R_GCPA_MF, 'gd-', 'Linewidth', lw);
% hold on;
% 
% grid on;
% legend('随机导频分配', '随机导频分配soft','加权图染色导频分配','加权图染色导频分配soft','图染色导频分配',  'Location', 'Best');
% xlabel('基站天线数 (M)');
% ylabel('用户下行平均可达速率 (Mbits/s/user)');
% hold off;

toc;















