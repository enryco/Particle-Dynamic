function [ vektorneu ] = koordinatenTrans( vektor, phi )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

vektorneu = [(vektor(1)*cos(phi) + vektor(2)*sin(phi))...
    (-vektor(1)*sin(phi) + vektor(2)*cos(phi))];



end

