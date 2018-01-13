clear;
clc; 
tic;

%% system parameters
M = 512;   
Gamma0 = 0.2;
% number of BS antennas M
L = 19;
K = 10;
S = 15;
lambda=0.1;

mu2=0.1;
J=floor(L*K*mu2);
%Test_num = 300;
Test_num = 15;
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



UL_SINR_CS_MF = zeros(K, L, Test_num);
UL_SINR_CSsoft_MF = zeros(K, L, Test_num);
UL_SINR_WGCPA_MF = zeros(K, L, Test_num);
UL_SINR_WGCPAsoft_MF = zeros(K, L, Test_num);
UL_SINR_WGCPAIM_MF = zeros(K, L, Test_num);


% DL_SINR_CS_MF = zeros(K, L, Test_num);
% DL_SINR_CSsoft_MF = zeros(K, L, Test_num);
% DL_SINR_WGCPA_MF = zeros(K, L, Test_num);
% DL_SINR_WGCPAsoft_MF = zeros(K, L, Test_num);
% DL_SINR_WGCPAIM_MF = zeros(K, L, Test_num);
% DL_SINR_GCPA_MF = zeros(K, L, Test_num);

% UL_R_CS_MF = zeros(K, L, Test_num);
% UL_R_WGCPA_MF = zeros(K, L, Test_num);
% UL_R_WGCPAIM_MF = zeros(K, L, Test_num);
% UL_R_GCPA_MF = zeros(K, L, Test_num);
% 
% 
% DL_R_CS_MF = zeros(K, L, Test_num);
% DL_R_WGCPA_MF = zeros(K, L, Test_num);
% DL_R_WGCPAIM_MF = zeros(K, L, Test_num);
% DL_R_GCPA_MF = zeros(K, L, Test_num);


%% Simulation
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
        
        UL_SINR_CS_MF(:,:,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);
               
    %    DL_SINR_CS_MF(:,:,i_test) = F_DLSINR_CAL(L, K, M, H, P, rho_pilot, rho_dl, 2);
        
      %% Random Pilot Assignment For Certen_SOFT
      
        [P eu] = F_CS_SOFT(L, K, S, Beta,lambda);
%         for i = 1:L                                                         % random pilot assignment
%             temp = randperm(S);
%             P(i,:) = temp(1:K);
%         end
        
        UL_SINR_CSsoft_MF(:,:,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);
               
      %  DL_SINR_CSsoft_MF(:,:,i_test) = F_DLSINR_CAL(L, K, M, H, P, rho_pilot, rho_dl, 2);

        
       %% Weighted Graph Coloring Based Pilot Assignment WGC-PA
        P = F_WGCPA_Pilot(L, K, S, Beta);
        
        UL_SINR_WGCPA_MF(:,:,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);

     %   DL_SINR_WGCPA_MF(:,:,i_test) = F_DLSINR_CAL(L, K, M, H, P, rho_pilot, rho_dl, 2);
       
       %% Weighted Graph Coloring Based Pilot Assignment For Certen_SOFT
        [P eu] = F_WGCPA_soft_Pilot(L, K, S, Beta,lambda);
        
        UL_SINR_WGCPAsoft_MF(:,:,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);

    %    DL_SINR_WGCPAsoft_MF(:,:,i_test) = F_DLSINR_CAL(L, K, M, H, P, rho_pilot, rho_dl, 2);
        
        
      %% IM_Weighted Graph Coloring Based Pilot Assignment WGC-PA
        P2 = F_WGCPA_Pilot_IM(L, K, S, Beta, J);
        
        UL_SINR_WGCPAIM_MF(:,:,i_test) = F_ULSINR_CAL(L, K, M, H, P2, rho_pilot, rho_ul, 1);



end

mu0=0.05;
T_1=1-S/K*mu0;
S_soft=S+eu*7;
T_soft=1-S_soft/K*mu0;
S_WIM=S+J;
T_WIM=1-S_WIM/K*mu0;

for i_test = 1:Test_num
    for l=1:L
        for k=1:K
     %       UL_SINR_CS_MF(k,l,i_test)=10^(UL_SINR_CS_MF(k,l,i_test)/10);
            UL_SINR_CS_MF(k,l,i_test)=Bw*T_1*log2(1+UL_SINR_CS_MF(k,l,i_test));
            %UL_SINR_CSsoft_MF(k,l,i_test)=Bw*T_ul/T_0*log2(1+UL_SINR_CSsoft_MF(k,l,i_test));
    %        UL_SINR_CSsoft_MF(k,l,i_test)=10^(UL_SINR_CSsoft_MF(k,l,i_test)/10);
            UL_SINR_CSsoft_MF(k,l,i_test)=Bw*T_soft*log2(1+UL_SINR_CSsoft_MF(k,l,i_test));
   %         DL_SINR_CS_MF(k,l,i_test)=10^(DL_SINR_CS_MF(k,l,i_test)/10);
   %         DL_SINR_CS_MF(k,l,i_test)=Bw*T_ul/T_0*log2(1+DL_SINR_CS_MF(k,l,i_test));
            %DL_SINR_CSsoft_MF(k,l,i_test)=Bw*T_ul/T_0*log2(1+DL_SINR_CSsoft_MF(k,l,i_test));
     %       DL_SINR_CSsoft_MF(k,l,i_test)=10^(DL_SINR_CSsoft_MF(k,l,i_test)/10);
     %       DL_SINR_CSsoft_MF(k,l,i_test)=Bw*T_soft*log2(1+DL_SINR_CSsoft_MF(k,l,i_test));
    %        UL_SINR_WGCPA_MF(k,l,i_test)=10^(UL_SINR_WGCPA_MF(k,l,i_test)/10);
            UL_SINR_WGCPA_MF(k,l,i_test)=Bw*T_1*log2(1+UL_SINR_WGCPA_MF(k,l,i_test));
            UL_SINR_WGCPAIM_MF(k,l,i_test)=Bw*T_WIM*log2(1+UL_SINR_WGCPAIM_MF(k,l,i_test));
           % UL_SINR_WGCPAsoft_MF(k,l,i_test)=Bw*T_ul/T_0*log2(1+UL_SINR_WGCPAsoft_MF(k,l,i_test));
    %       UL_SINR_WGCPAsoft_MF(k,l,i_test)=10^(UL_SINR_WGCPAsoft_MF(k,l,i_test)/10);
            UL_SINR_WGCPAsoft_MF(k,l,i_test)=Bw*T_soft*log2(1+UL_SINR_WGCPAsoft_MF(k,l,i_test));
    %        DL_SINR_WGCPA_MF(k,l,i_test)=10^(DL_SINR_WGCPA_MF(k,l,i_test)/10);
    %        DL_SINR_WGCPA_MF(k,l,i_test)=Bw*T_ul/T_0*log2(1+ DL_SINR_WGCPA_MF(k,l,i_test));
           % DL_SINR_WGCPAsoft_MF(k,l,i_test)=Bw*T_ul/T_0*log2(1+ DL_SINR_WGCPAsoft_MF(k,l,i_test));
    %       DL_SINR_WGCPAsoft_MF(k,l,i_test)=10^(DL_SINR_WGCPAsoft_MF(k,l,i_test)/10);
   %        DL_SINR_WGCPAsoft_MF(k,l,i_test)=Bw*T_soft*log2(1+ DL_SINR_WGCPAsoft_MF(k,l,i_test));
   %        UL_SINR_GCPA_MF(k,l,i_test)=10^(UL_SINR_GCPA_MF(k,l,i_test)/10);
  %         DL_SINR_GCPA_MF(k,l,i_test)=10^(DL_SINR_GCPA_MF(k,l,i_test)/10);
   %         DL_SINR_GCPA_MF(k,l,i_test)=Bw*T_ul/T_0*log2(1+DL_SINR_GCPA_MF(k,l,i_test));
        end
    end
end
Tau = 20;

%% UL SINR
figure; 

y = reshape(UL_SINR_CS_MF,1,L*K*Test_num);
ymin=min(y); 
ymax=max(y); 
x=linspace(ymin,ymax,Tau); 
yy=hist(y,x);
yy=yy/length(y)/(x(2)-x(1));
s=0; 
for i=2:length(x) 
    s=[s,trapz(x([1:i]),yy([1:i]))]; 
end 
plot(x,s,'bo-');
hold on;

y = reshape(UL_SINR_CSsoft_MF,1,L*K*Test_num);
ymin=min(y); 
ymax=max(y); 
x=linspace(ymin,ymax,Tau); 
yy=hist(y,x);
yy=yy/length(y)/(x(2)-x(1));
s=0; 
for i=2:length(x) 
    s=[s,trapz(x([1:i]),yy([1:i]))]; 
end 
plot(x,s,'k>-');
hold on;


y = reshape(UL_SINR_WGCPA_MF,1,L*K*Test_num);
ymin=min(y); 
ymax=max(y); 
x=linspace(ymin,ymax,Tau); 
yy=hist(y,x);
yy=yy/length(y)/(x(2)-x(1));
s=0; 
for i=2:length(x) 
    s=[s,trapz(x([1:i]),yy([1:i]))]; 
end 
plot(x,s,'r^-');
hold on;
%{
y = reshape(UL_SINR_WGCPAIM_MF,1,L*K*Test_num);
ymin=min(y); 
ymax=max(y); 
x=linspace(ymin,ymax,Tau); 
yy=hist(y,x);
yy=yy/length(y)/(x(2)-x(1));
s=0; 
for i=2:length(x) 
    s=[s,trapz(x([1:i]),yy([1:i]))]; 
end 
plot(x,s,'mx-');
hold on;
%}
y = reshape(UL_SINR_WGCPAsoft_MF,1,L*K*Test_num);
ymin=min(y); 
ymax=max(y); 
x=linspace(ymin,ymax,Tau); 
yy=hist(y,x);
yy=yy/length(y)/(x(2)-x(1));
s=0; 
for i=2:length(x) 
    s=[s,trapz(x([1:i]),yy([1:i]))]; 
end 
plot(x,s,'cs-');
hold on;

% y = 10*log10(reshape(UL_SINR_WGCPAIM_MF,1,L*K*Test_num));
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'ys-');
% hold on;




grid on;
%axis([-80 40 0 1]);
legend('Random','SPRS','WGC-PD','SPRS+WGC-PD', 'Location', 'Best');
xlabel('Average uplink achievable rate per user C (bps/Hz)');
ylabel('CDF of the users uplink achievable rate');
hold off;


%% DL SINR
% figure; 
% 
% y = reshape(DL_SINR_CS_MF,1,L*K*Test_num);
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'bo-');
% hold on;
% 
% y = reshape(DL_SINR_CSsoft_MF,1,L*K*Test_num);
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'k>-');
% hold on;
% 
% y = reshape(DL_SINR_WGCPA_MF,1,L*K*Test_num);
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'r^-');
% hold on;
% 
% y = reshape(DL_SINR_WGCPAsoft_MF,1,L*K*Test_num);
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'cs-');
% hold on;
% 
% % y = 10*log10(reshape(DL_SINR_WGCPAIM_MF,1,L*K*Test_num));
% % ymin=min(y); 
% % ymax=max(y); 
% % x=linspace(ymin,ymax,Tau); 
% % yy=hist(y,x);
% % yy=yy/length(y)/(x(2)-x(1));
% % s=0; 
% % for i=2:length(x) 
% %     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% % end 
% % plot(x,s,'ys-');
% % hold on;
% 
% 
% y = reshape(DL_SINR_GCPA_MF,1,L*K*Test_num);
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'gd-');
% hold on;
% 
% grid on;
% %axis([-80 40 0 1]);
% legend('随机导频分配','随机导频分配soft', '加权图染色导频分配','加权图染色导频分配soft', '图染色导频分配', 'Location', 'Best');
% xlabel('用户下行平均可达速率');
% ylabel('累积分布函数 (CDF)');
% hold off;



toc;
%}














