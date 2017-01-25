clear all
close all
clc

lambdaratio = 2;
deltaratio = 2;

OLD = [0.460 7.096 57.720 450.225 3400.342];
    
NEW = [0.680 3.257 9.580 32.244 141.227];

INTERCELLS = [195 768 1534 3073 6161];

subplot(1,2,1)
WORSTCASE = [0.680 .680*(INTERCELLS(2)/INTERCELLS(1))*(lambdaratio^2+log(deltaratio^1)) .680*(INTERCELLS(3)/INTERCELLS(1))*((lambdaratio^2)^2+log(deltaratio^2)) .680*(INTERCELLS(4)/INTERCELLS(1))*((lambdaratio^3)^2+log(deltaratio^3)) .680*(INTERCELLS(5)/INTERCELLS(1))*((lambdaratio^4)^2+log(deltaratio^4))];
loglog(INTERCELLS,OLD,'-or');
title('One Intersection Loop')
xlabel('Number of Intersecting Edges')
ylabel('Seconds')
p1 = polyfit(log(INTERCELLS),log(OLD),1)
hold on
loglog(INTERCELLS,NEW, '--s')
p2 = polyfit(log(INTERCELLS),log(NEW),1)
p3 = polyfit(log(INTERCELLS),log(WORSTCASE),1)
%loglog(INTERCELLS,WORSTCASE, '-*','color','black')
legend('VTK 6.2.0 Implementation','New Implementation','Worst Case')
str1 = sprintf('    m = %.4f',p1(1));
str2 = sprintf('    m = %.4f',p2(1));
str3 = sprintf('m = %.4f   ',p3(1));
text(INTERCELLS(4),OLD(4),str1,'color','red')
text(INTERCELLS(2),NEW(2),str2,'color','blue')
%text(INTERCELLS(3),WORSTCASE(3),str3,'HorizontalAlignment','right')

lambdaratio = 2;
deltaratio = 4;

OLD = [1.194 2.335 6.733 38.934 284.667];
    
NEW = [2.314 3.516 6.979 21.688 69.285];

INTERCELLS = [1179 1869 3244 6738 13272];

WORSTCASE = [2.314 2.314*(INTERCELLS(2)/INTERCELLS(1))*(lambdaratio^2+log(deltaratio^1)) 2.314*(INTERCELLS(3)/INTERCELLS(1))*((lambdaratio^2)^2+log(deltaratio^2)) 2.314*(INTERCELLS(4)/INTERCELLS(1))*((lambdaratio^3)^2+log(deltaratio^3)) 2.314*(INTERCELLS(5)/INTERCELLS(1))*((lambdaratio^4)^2+log(deltaratio^4))];

subplot(1,2,2)
loglog(INTERCELLS,OLD,'-or');
title('Four Intersection Loops')
xlabel('Number of Intersecting Edges')
ylabel('Seconds')
axis([10^3 10^4.2 1 10^3])
p1 = polyfit(log(INTERCELLS),log(OLD),1)
hold on
loglog(INTERCELLS,NEW, '--s')
p2 = polyfit(log(INTERCELLS),log(NEW),1)
p3 = polyfit(log(INTERCELLS),log(WORSTCASE),1)
%loglog(INTERCELLS,WORSTCASE, '-*','color','black')
legend('VTK 6.2.0 Implementation','New Implementation','Worst Case')
str1 = sprintf('m = %.4f   ',p1(1));
str2 = sprintf('    m = %.4f',p2(1));
str3 = sprintf('m = %.4f   ',p3(1));
text(INTERCELLS(4),OLD(4),str1,'color','red','HorizontalAlignment','right')
text(INTERCELLS(4),NEW(4),str2,'color','blue')
%text(INTERCELLS(3),WORSTCASE(3),str3,'HorizontalAlignment','right')

