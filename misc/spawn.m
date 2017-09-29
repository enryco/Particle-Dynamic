%SPawn Script

N = guidata(numberOfObjectsEdit);
% N = get(numberOfObjectsEdit, 'String');
N = str2double(N);

atoms = atomClass;

hold on

for i=1:N
    atoms(i) = atomClass;
    atoms(i).coordinates = [rand rand];
    atoms(i).radius = rand;
   
    viscircles([atoms(i).coordinates(1)...
        atoms(i).coordinates(2)],...
        A(i).radius);
end