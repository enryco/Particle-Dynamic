%% Initialization
clear all;
close all;
clc;
clf;
cla;

%% Box Size and Axis Initialization

boxSize = [1000 1000];   %[width hight] of the Box

axis equal;
axis([0 boxSize(1) 0 boxSize(2)]);
hold on

%% Declaring Variables

numberOfAtoms = 10;             %Number of Atoms
radius = 71;      %Radius of Atoms
itMax = 10000;      %Maximum iterations for spawning atoms
dt = 0.01;          %Time
velocity = 130;     %Velocity
acceleration = 0;   %Acceleration - will be calculated later via F = m*a

%% Initialize First Atom as Object
atoms = atomClass;

%% Spawn Atoms
%Create and randomly Spawn other AtomObjects randoy


for i=1:numberOfAtoms
    %% Create Atom Objects & Set atributs
    
    
    
    atoms(i) = atomClass;       %Create object of Class 'atomClass'
    atoms(i).radius = radius;   %assign radii
    atoms(i).direction = [rand-0.5 rand-0.5];   %randomly assign direction of atom movement
    atoms(i).direction = atoms(i).direction/norm(atoms(i).direction); %normalize direction vector
    atoms(i).coordinates = [rand*boxSize(1) rand*boxSize(2)];
    
    atoms(i).velocity = velocity*atoms(i).direction;
    atoms(i).acceleration = acceleration*atoms(i).direction;
    
    %% Check if Atom is placed NOT inside another atom:
    checkiterations = 0; %Initializing to count iterations
    if i > 1
        check = false;
    else
        check = true;
    end
    while ~check
        checkiterations = checkiterations+1; %counting iterations to break
        
        % Ckeck wether Atoms overlap or are outside the box:
        for j=1:i-1
            if atoms(1).overlap(atoms([i j])) %check overlap
                %  || atoms(1).isInBox(atoms(i),boxSize) == false
                atoms(i).coordinates = [rand*boxSize(1) rand*boxSize(2)]; %set new coordinates
                %disp(['new corrdinates ' num2str(atoms(i).coordinates)])
                break
            end
            
            %Check if last atom was checked
            if j == i-1
                check = true; %set check to true as all atoms were checked with no overlapping
                %'done'
            end
        end
        
        %Break conditions
        if checkiterations > itMax
            atoms(i) = [];
            disp('There might not be enough space to place all atoms')
            break
        end
    end
    
    %% break conditions
    if checkiterations > itMax
        disp(['Spawn stopped at ' num2str(i-1) ' Atoms'])
        break
    end
    
    %% give atom a color based on position
    if atoms(i).coordinates(1) < boxSize(1)/2
        atoms(i).color = 'red';
    else
        atoms(i).color = 'blue';
    end
    
    %% Plotting
%     plot(atoms(i).coordinates(1),...
%         atoms(i).coordinates(2),'o','Color',atoms(i).color);
            viscircles([atoms(i).coordinates(1)...
                atoms(i).coordinates(2)],...
                atoms(i).radius);
end

pause

%% Let The Atoms Fly!

%{
%initialize first two iterations
% X(1) = x0;
% X(2) = X(1) + v0*(tStart+dt) + 0.5*a*(tStart+dt)^2;

for i=1:size(atoms,2)
    atoms(i).coordinates = [atoms(i).coordinates;...
        atoms(i).coordinates + velocityStart*dt*atoms(i).direction + 0.5*acceleration*dt^2];
end
%}

iterations = 0;
accelerationArrayX = zeros(size(atoms,2));
accelerationArrayY = zeros(size(atoms,2));
velocityDiff = 0;
velocitySumStart = [0 0];

for i=1:size(atoms,2)
    velocitySumStart = velocitySumStart + atoms(i).velocity;
end

%dn(velocitySumStart)
%dn(norm(velocitySumStart))

% close all

while 1
    diffusionRightSide = 0;
    diffusionLeftSide = 0;
    
    velocitySum = 0;
    velocityDiff = 0;
    
    iterations = iterations + 1;
    
    cla
    
    %Calc new accelerations
    for i=1:(size(atoms,2)-1)
        for j=(i+1):size(atoms,2)
            [accAtom1,accAtom2] =  aLJP(atoms(i),atoms(j));
            accelerationArrayX(i,j) = accAtom1(1);
            accelerationArrayX(j,i) = accAtom2(1);
            accelerationArrayY(i,j) = accAtom1(2);
            accelerationArrayY(j,i) = accAtom2(2);
        end
    end
    
    for i=1:size(atoms,2)
        atoms(i).acceleration = [sum(accelerationArrayX(i,:)) sum(accelerationArrayY(i,:))];
    end
    %
    %     disp('AccCalc Completed')
    
    
    %Updating Positions
    for i=1:size(atoms,2)
        atoms(i).coordinates = atoms(i).coordinates +... %set new coordinates
            atoms(i).velocity*dt +...
            0.5*atoms(i).acceleration*(dt^2);
        %velocity = velocity + acceleration*dt; %set new velocity
        atoms(i).velocity = atoms(i).velocity + atoms(i).acceleration*dt;
        
        
        %put atom inside box if its outside
        if atoms(i).canPutBackIntoBox(atoms(i),boxSize) %check wether it's possible to put atom back inside box
            atoms(i).coordinates = atoms(i).putBackIntoBox(atoms(i),boxSize);
        else
            disp(['Iterations' num2str(iterations)])
            error(['Atom ' num2str(i) ' velocities exeed box marging' 'At:' num2str(atoms(i).coordinates)])
            
        end
        
        %count atoms on the left side
        if strcmp(atoms(i).color,'red') && atoms(i).coordinates(1) < boxSize(1)/2
            diffusionLeftSide = diffusionLeftSide + 1;
        elseif strcmp(atoms(i).color,'red') && atoms(i).coordinates(1) > boxSize(1)/2
            diffusionRightSide = diffusionRightSide + 1;
        end
        
%                 plot(atoms(i).coordinates(1),...
%                     atoms(i).coordinates(2),'o','Color',atoms(i).color);
                        viscircles([atoms(i).coordinates(1)...
                            atoms(i).coordinates(2)],...
                            atoms(i).radius);
        % function to check if it is inside box
        
        
        
        velocitySum = velocitySum + atoms(i).velocity;
        
    end
    % dn(velocitySum)
    velocityDiff = (velocitySumStart - velocitySum)/size(atoms,2);
    % dn(velocityDiff)
    
    for i=1:size(atoms,2)
        atoms(i).velocity = atoms(i).velocity + velocityDiff;
    end
    
    velocitySum = 0;
    for i=1:size(atoms,2)
        velocitySum = velocitySum + atoms(i).velocity;
    end
    %dn(norm(velocitySum))
    % plot(iterations,diffusionLeftSide)
    %         disp(['Diff L:' num2str(diffusionLeftSide)...
    %             'Diff R:' num2str(diffusionRightSide) 'sum ' num2str(diffusionLeftSide+diffusionRightSide)])
    
        pause(0.02)
    if mod(iterations,100) == 0
        disp(i)
    end
end
