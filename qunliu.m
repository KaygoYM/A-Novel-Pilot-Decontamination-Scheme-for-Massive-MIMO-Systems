function qunliu( xx,yy,r,ceng,bp,up )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hang=linspace(0,2*pi,7);
    for i=1:1:ceng
        for j=1:1:length(hang)  
        theta(i,i*j-i+1) = hang(j);
        cx(i,i*j-i+1)=sqrt(3)*i*r*cos(theta(i,i*j-i+1));  %每层上角度为pi/3,2*pi/3...小区中心
        cy(i,i*j-i+1)=sqrt(3)*i*r*sin(theta(i,i*j-i+1));
       
        end
    end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if(ceng>=2)
         for i=2:1:ceng
                for j=1:1:length(hang)-1                  
                        for k=1:1:i-1     %求那些非pi/3整数倍角度上的中心，依据第n层的是n分点
                       cx(i,i*j-i+1+k)=((i-k)*cx(i,i*j-i+1)+k*cx(i,i*j+1))/i;
                       cy(i,i*j-i+1+k)=((i-k)*cy(i,i*j-i+1)+k*cy(i,i*j+1))/i;
                        end
                end
         end
     end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zx=cx+xx;
zy=cy+yy;
for i=1:1:ceng
    for j=1:1:6*i
    cirliu(zx(i,j),zy(i,j),r,bp,up);
    hold on;
    end
end