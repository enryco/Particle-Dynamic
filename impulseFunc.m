function [v1neu, v2neu] = impulseFunc(...
    coordinatesAtom1,coordinatesAtom2,...
    velocityAtom1,velocityAtom2)
%Funktion zur Berechnung eines Stoﬂes

z = coordinatesAtom2 - coordinatesAtom1; %zentralvektor 
y = [0 1];
phi = acos( ( y*z' ) / ( sqrt(y*y')*sqrt(z*z') ) );

v1z = koordinatenTrans(velocityAtom1, phi);
v2z = koordinatenTrans(velocityAtom2, phi);

v1neuz = [v1z(1) v2z(2)];
v2neuz = [v2z(1) v1z(2)];

v1neu = koordinatenTrans(v1neuz, -phi);
v2neu = koordinatenTrans(v2neuz, -phi);

end