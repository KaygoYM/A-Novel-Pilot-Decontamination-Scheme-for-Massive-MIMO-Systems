function cirliu(xx,yy,r,bp,up)
theta = linspace(pi/6,13*pi/6,7);
plot(r*cos(theta)+xx,r*sin(theta)+yy,'g-','LineWidth',2);
hold on;

for i=1:19
    plot(bp(i,1),bp(i,2),'k^','MarkerFaceColor','k');
    hold on
end
 
for l=1:19
    for k=1:5
    plot(up(l,k,1),up(l,k,2),'rx');
    hold on
    end 
end

axis([-5500 5500 -5000 5000]);
%legend('小区边界', '基站','用户',  'Location', 'Best');
