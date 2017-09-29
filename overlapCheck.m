function overlapBool = overlapCheck(coordinatesAtom1,coordinatesAtom2,radiusAtom1,radiusAtom2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 %calc distance between atoms:
            distance = sqrt(sum((coordinatesAtom1 - coordinatesAtom2).^2));
            if distance <= (radiusAtom1 + radiusAtom2)
                overlapBool = true;
            else
                overlapBool = false;
            end

end

