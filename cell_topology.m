clear;
clc;

L=19; 
K=5;
R_Cell=1000;
r_Min = 100;
BaseP = zeros(L,2);                                                         % coords of BS
UserP = zeros(L,K,2); 
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
        r=R_Cell/cos(pi/6);
        
        for i=1:L
            for k=1:K
                x2 = random('unif', -R_Cell, R_Cell, 1, 1);
                y2 = random('unif', -R_Cell, R_Cell, 1, 1);
                while ((x2)^2+(y2)^2>R_Cell^2) || ((x2)^2+(y2)^2<r_Min^2)
                x2 = random('unif', -R_Cell, R_Cell, 1, 1);
                y2 = random('unif', -R_Cell, R_Cell, 1, 1);
                end
                x2 = x2 + BaseP(i,1);
                y2 = y2 + BaseP(i,2);
                UserP(i,k,1)=x2;
                UserP(i,k,2)=y2;
            end
        end

        qunliu(0,0,r,2,BaseP,UserP)               %第一个和第二个参数为最中心小区的横纵坐标，第三个参数为小区半径，第四个为所需画的层数

