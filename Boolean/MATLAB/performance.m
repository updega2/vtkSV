clear all
close all
clc

lambdaratio = 2;
deltaratio = 2;

OLD = [0.460 7.096 57.720 450.225 3400.342];
    
NEW = [0.680 3.257 9.580 32.244 141.227];

INTERCELLS = [195 768 1534 3073 6161];

WORSTCASE = [0.680 .680*(INTERCELLS(2)/INTERCELLS(1))*(lambdaratio^2+log(deltaratio^1)) .680*(INTERCELLS(3)/INTERCELLS(1))*((lambdaratio^2)^2+log(deltaratio^2)) .680*(INTERCELLS(4)/INTERCELLS(1))*((lambdaratio^3)^2+log(deltaratio^3)) .680*(INTERCELLS(5)/INTERCELLS(1))*((lambdaratio^4)^2+log(deltaratio^4))];
loglog(INTERCELLS,OLD,'-or');
title('Run Time for Boolean Operations')
xlabel('Number of Intersecting Edges')
ylabel('Seconds')
p1 = polyfit(log(INTERCELLS),log(OLD),1)
hold on
loglog(INTERCELLS,NEW, '--s')
p2 = polyfit(log(INTERCELLS),log(NEW),1)
p3 = polyfit(log(INTERCELLS),log(WORSTCASE),1)
loglog(INTERCELLS,WORSTCASE, '-*','color','black')
legend('Old Implementation','New Implementation','Worst Case')

