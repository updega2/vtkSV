clear all
close all
clc

lambdaratio = 2;
deltaratio = 4;

OLD = [1.194 2.335 6.733 38.934 284.667];
    
NEW = [2.314 3.516 6.979 21.688 69.285];

INTERCELLS = [1179 1869 3244 6738 13272];

WORSTCASE = [2.314 2.314*(INTERCELLS(2)/INTERCELLS(1))*(lambdaratio^2+log(deltaratio^1)) 2.314*(INTERCELLS(3)/INTERCELLS(1))*((lambdaratio^2)^2+log(deltaratio^2)) 2.314*(INTERCELLS(4)/INTERCELLS(1))*((lambdaratio^3)^2+log(deltaratio^3)) 2.314*(INTERCELLS(5)/INTERCELLS(1))*((lambdaratio^4)^2+log(deltaratio^4))];

loglog(INTERCELLS,OLD,'-or');
title('Run Time for Boolean Operations')
xlabel('Number of Intersecting Edges')
ylabel('Seconds')
axis([10^3 10^4.2 1 10^3])
p1 = polyfit(log(INTERCELLS),log(OLD),1)
hold on
loglog(INTERCELLS,NEW, '--s')
p2 = polyfit(log(INTERCELLS),log(NEW),1)
p3 = polyfit(log(INTERCELLS),log(WORSTCASE),1)
loglog(INTERCELLS,WORSTCASE, '-*','color','black')
legend('Old Implementation','New Implementation','Worst Case')

