function [ vektorneu ] = koordinatenTrans( vektor, phi )
%Transformieren der Koordinaten

vektorneu = [(vektor(1)*cos(phi) + vektor(2)*sin(phi))...
    (-vektor(1)*sin(phi) + vektor(2)*cos(phi))];

end

