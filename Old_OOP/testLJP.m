%% Test Script
clear all;
close all;
clc;
clf;
cla;

% axis manual
axis equal
% axis([-70 140 0 140]);
dt = 1;

atoms = atomClass;
atoms(2) = atomClass;

velocity = 1;
atoms(1).coordinates = [0 0];
atoms(2).coordinates = [340/6 10];
atoms(1).direction = [1 1];
atoms(2).direction = [-1 0.5];
atoms(1).velocity = velocity*atoms(1).direction;
atoms(2).velocity = velocity*atoms(2).direction;
atoms(1).acceleration = [0 0];
atoms(2).acceleration = [0 0];
atoms(1).color = 'red';
atoms(2).color = 'blue';
atoms(1).radius = 71;
atoms(2).radius = 71;

hold on
for i=1:2
    plot(atoms(i).coordinates(1),...
        atoms(i).coordinates(2),'o');
end


pause
% axis manual

% Gesamtenergie muss konstant bleiben...
iterations = 0;
while 1
    %     cla
    iterations = iterations+1;
    [atoms(1).acceleration,atoms(2).acceleration] = aLJP(atoms(1),atoms(2));
    %     atoms(1).acceleration
    %     atoms(2).acceleration
    
    for i=1:2
        atoms(i).coordinates = atoms(i).coordinates +... %set new coordinates
            atoms(i).velocity*dt +...
            0.5*atoms(i).acceleration*(dt^2);
        %velocity = velocity + acceleration*dt; %set new velocity
        atoms(i).velocity = atoms(i).velocity + atoms(i).acceleration*dt;
        plot(atoms(i).coordinates(1),...
            atoms(i).coordinates(2),'.','Color',atoms(i).color);
        %         viscircles([atoms(i).coordinates(1)...
        %             atoms(i).coordinates(2)],...
        %             atoms(i).radius);
    end
    
    if mod(iterations,1000) == 0
        pause(0.02)
    end
end
