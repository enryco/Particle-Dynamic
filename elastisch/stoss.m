function [v1neu, v2neu] = stoss(p1,p2,v1,v2)
%Funktion zur Berechnung eines Stoﬂes

z = p2 - p1; %zentralvektor 
y = [0 1];
phi = acos( ( y*z' ) / ( sqrt(y*y')*sqrt(z*z') ) );

v1z = koordinatenTrans(v1, phi);
v2z = koordinatenTrans(v2, phi);

v1neuz = [v1z(1) v2z(2)];
v2neuz = [v2z(1) v1z(2)];

v1neu = koordinatenTrans(v1neuz, -phi);
v2neu = koordinatenTrans(v2neuz, -phi);

end