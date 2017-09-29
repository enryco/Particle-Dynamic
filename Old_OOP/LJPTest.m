%% Test of The Leonard Jones Potential

sigma = 340.5; %pm
epsilon = 1.653e-21; % 10-21 J

x1 = 0;
y1 = 0;
x2 = 1.122*340
y2 = 0;%  1e-10;


Fx1 = 	-4*epsilon*((12*sigma^12*(x2-x1))/((x2-x1)^2+(y2-y1)^2)^7-(6*sigma^6*(x2-x1))/((x2-x1)^2+(y2-y1)^2)^4);
Fy1 =	-4*epsilon*((12*sigma^12*(y2-y1))/((x2-x1)^2+(y2-y1)^2)^7-(6*sigma^6*(y2-y1))/((x2-x1)^2+(y2-y1)^2)^4);
Fx2 =	-4*epsilon*((6*sigma^6*(x2-x1))/((x2-x1)^2+(y2-y1)^2)^4-(12*sigma^12*(x2-x1))/((x2-x1)^2+(y2-y1)^2)^7);
Fy2 =	-4*epsilon*((6*sigma^6*(y2-y1))/((x2-x1)^2+(y2-y1)^2)^4-(12*sigma^12*(y2-y1))/((x2-x1)^2+(y2-y1)^2)^7);


[Fx1 Fy1]
[Fx2 Fy2]

hold on
plot(x1,y1,'+',x2,y2,'+')
line([x1 Fx1],[y1 Fy1])
