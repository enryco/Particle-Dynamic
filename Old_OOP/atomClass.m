classdef atomClass
    %atomClass is a class representing an atom
    %   Detailed explanation goes here
    
    properties
        name
        radius
        mass
        velocity
        acceleration
        coordinates
        direction
        color
    end
    
    methods (Static)
        
        %Function to determin wether two atoms overlap or not
        %outputs a bool
        function overlapBool = overlap(twoAtomsArray)
            %calc distance between atoms:
            distance = sqrt(sum((twoAtomsArray(1).coordinates - twoAtomsArray(2).coordinates).^2));
            if distance < (twoAtomsArray(1).radius + twoAtomsArray(2).radius)
                overlapBool = true;
            else
                overlapBool = false;
            end
        end
        
        %function to check wether atom is still inside box
        %true if is in box; fals if it is Not in box
        function isInBoxBool = isInBox(atom, boxSize)
            if (atom.coordinates(1) - atom.radius) < 0 ... %check left margin
                    || (atom.coordinates(1) + atom.radius) > boxSize(1) ... %check right margin
                    || (atom.coordinates(2) - atom.radius) < 0 ... %check bottom margin
                    || (atom.coordinates(2) + atom.radius) > boxSize(2) %check top margin
                isInBoxBool = false;
            else
                isInBoxBool = true;
                %disp('outside box')
            end
            
        end
        
       %{
function newCoordinates = putBackIntoBox(atom, boxSize)
            new_x = 0;
            new_y = 0;
            
            if atom.coordinates(1) - atom.radius < 0 %check left margin
                new_x = atom.coordinates(1)+boxSize(1);
                newCoordinates = [new_x atom.coordinates(2)];
                
            elseif atom.coordinates(1) + atom.radius > boxSize(1) %check right margin
                new_x = atom.coordinates(1) - boxSize(1);
                newCoordinates = [new_x atom.coordinates(2)];
                
            elseif atom.coordinates(2) - atom.radius < 0 %check bottom margin
                new_y = atom.coordinates(2)+boxSize(2);
                newCoordinates = [atom.coordinates(1) new_y];
                
            elseif atom.coordinates(2) + atom.radius > boxSize(2) %check top margin
                new_y = atom.coordinates(2)-boxSize(2);
                newCoordinates = [atom.coordinates(1) new_y];
            else
                newCoordinates = atom.coordinates;
            end
            
        end
        %}
        
        function newCoordinates = putBackIntoBox(atom, boxSize)
           
            newCoordinates = mod(atom.coordinates,boxSize);
            
        end
        
         function canPutBackIntoBoxBool = canPutBackIntoBox(atom, boxSize)
            if (atom.coordinates(1) - atom.radius) < 0-boxSize(1) ... %check left margin
                    || (atom.coordinates(1) + atom.radius) > 2*boxSize(1) ... %check right margin
                    || (atom.coordinates(2) - atom.radius) < 0-boxSize(2) ... %check bottom margin
                    || (atom.coordinates(2) + atom.radius) > 2*boxSize(2) %check top margin
                canPutBackIntoBoxBool = false;
            else
                canPutBackIntoBoxBool = true;
                %                 disp('outside box')
            end
            
         end
        
    end
end

