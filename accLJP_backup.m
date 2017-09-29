function [accAtom1,accAtom2] = accLJP(coordinatesAtom1,coordinatesAtom2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

sigma = 340.5; %pm
epsilon = 1.653e-21; % 10-21 J
mass = 6.634e-26*1e12; %kg -> pg  (kilograms) (39.948u)

x1 = coordinatesAtom1(1);
y1 = coordinatesAtom1(2);
x2 = coordinatesAtom2(1);
y2 = coordinatesAtom2(2);

Fx1 = 	-4*epsilon*((12*sigma^12*(x2-x1))/((x2-x1)^2+(y2-y1)^2)^7-(6*sigma^6*(x2-x1))/((x2-x1)^2+(y2-y1)^2)^4);
Fy1 =	-4*epsilon*((12*sigma^12*(y2-y1))/((x2-x1)^2+(y2-y1)^2)^7-(6*sigma^6*(y2-y1))/((x2-x1)^2+(y2-y1)^2)^4);
Fx2 =	-4*epsilon*((6*sigma^6*(x2-x1))/((x2-x1)^2+(y2-y1)^2)^4-(12*sigma^12*(x2-x1))/((x2-x1)^2+(y2-y1)^2)^7);
Fy2 =	-4*epsilon*((6*sigma^6*(y2-y1))/((x2-x1)^2+(y2-y1)^2)^4-(12*sigma^12*(y2-y1))/((x2-x1)^2+(y2-y1)^2)^7);

accAtom1 = [Fx1 Fy1]/mass;
accAtom2 = [Fx2 Fy2]/mass;

end

