function [v1neu, v2neu] = impulseFunc(...
    coordinatesAtom1,coordinatesAtom2,...
    velocityAtom1,velocityAtom2)
%Funktion zur Berechnung eines Stoßes

z = coordinatesAtom2 - coordinatesAtom1; %zentralvektor 
y = [0 1];  %y-achsen vektor
phi = acos( ( y*z' ) / ( sqrt(y*y')*sqrt(z*z') ) ); %calc phi

%Koordinaten Transfer
v1z = koordinatenTrans(velocityAtom1, phi); 
v2z = koordinatenTrans(velocityAtom2, phi); 

%Berechnung von v1neu und v2neu
v1neuz = [v1z(1) v2z(2)];  
v2neuz = [v2z(1) v1z(2)];

%Zurückrechnen ins eigentlich Koordinaten-System
v1neu = koordinatenTrans(v1neuz, -phi);
v2neu = koordinatenTrans(v2neuz, -phi);

end