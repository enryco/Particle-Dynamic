function [accAtom1,accAtom2] = accLJP(coordinatesAtom1,coordinatesAtom2,sigma6,sigma12,epsilon,mass)
%Funktion zum Berechnen der Kraft des LJP

x1 = coordinatesAtom1(1);
y1 = coordinatesAtom1(2);
x2 = coordinatesAtom2(1);
y2 = coordinatesAtom2(2);

deltax = x2-x1;
deltax2 = deltax^2;
deltay = y2-y1;
deltay2 = deltay^2;

Fx1 = 	-4*epsilon*((12*sigma12*(deltax))/(deltax2+deltay2)^7-(6*sigma6*(deltax))/(deltax2+deltay2)^4);
Fy1 =	-4*epsilon*((12*sigma12*(deltay))/(deltax2+deltay2)^7-(6*sigma6*(deltay))/(deltax2+deltay2)^4);
Fx2 =	-Fx1;
Fy2 =	-Fy1;

accAtom1 = [Fx1 Fy1]/mass;
accAtom2 = [Fx2 Fy2]/mass;

end

