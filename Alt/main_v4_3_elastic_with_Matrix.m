%% Initialization
clear all;
close all;
% clc;
clf;
cla;

%% Box Size and Axis Initialization

boxSize = 1e4 * [1 1];   %[width hight] of the Box

axis equal;
axis([0 boxSize(1) 0 boxSize(2)]);
hold on

%% Declaring Variables

numberOfAtoms = 20; %Number of Atoms
radius = 71;        %Radius of Atoms
itMax = 10000;      %Maximum iterations for spawning atoms
dt = 0.01;             %Time
velocity = 416*1e12;      %Velocity
acceleration = 0;   %Acceleration - will be calculated later via F = m*a

dt = radius/velocity

%% jetzt gehts los

radii = ones(numberOfAtoms,1)*radius;

%% Spawning inside box without overlaping box margins
% coordinates = [(rand(numberOfAtoms,1)*(boxSize(1)-2*radius))+radius...
%     (rand(numberOfAtoms,1)*(boxSize(2)-2*radius))+radius];

coordinates = zeros(numberOfAtoms,2);   %precreate coordinates matrix
coordinates(1,:) = [(rand*(boxSize(1)-2*radius))+radius... %asign first atom
    (rand*(boxSize(2)-2*radius))+radius];

%spawn other atoms while checking if newly spawned atom overlaps with any
%otehr atom
for i=2:numberOfAtoms
    
    %create coordinates:
    coordinates(i,:) = [(rand*(boxSize(1)-2*radius))+radius...
        (rand*(boxSize(2)-2*radius))+radius];
    
    %check wether newly created atom overlaps:
    check = false;
    checkIterations = 0;
    while ~check
        checkIterations = checkIterations + 1;
        for j=1:i-1
            if overlapCheck(coordinates(i,:),coordinates(j,:),radii(i),radii(j))
                coordinates(i,:) = [(rand*(boxSize(1)-2*radius))+radius...
                    (rand*(boxSize(2)-2*radius))+radius];
                break
            end
            if j==i-1
                check = true;
            end
        end
        if checkIterations > 10000
            error('Spawning Atoms without overlap might not possible')
        end
    end
    
end

leftSide = coordinates(:,1) > boxSize(1)/2;
circleColor = [leftSide zeros(numberOfAtoms,1) not(leftSide)];


%% setting random directions and apply velocity
direction = rand(numberOfAtoms,2)-0.5; %assign random directions
directionNorm = sqrt(direction(:,1).^2 + direction(:,2).^2); %calc normVector
direction = direction./[directionNorm directionNorm]; %normalizing direction
velocities = velocity*direction; %assign velocity


leftSideIndex = find(leftSide==1);
rightSideIndex = find(leftSide==0);
viscircles(coordinates(leftSideIndex,:),radii(leftSideIndex),'EdgeColor','b');
viscircles(coordinates(rightSideIndex,:),radii(rightSideIndex),'EdgeColor','r');

pause

%% let the atoms fly
t = 0;
while 1
    t = t+1;
    cla
    
    dt = radius/max(sqrt(velocities(:,1).^2+velocities(:,2).^2))
    
    
    %calculate and update new velocities
    
    X = coordinates(:,1)*ones(1,numberOfAtoms);
    X = triu(abs(X-X'));
    Y = coordinates(:,2)*ones(1,numberOfAtoms);
    Y = triu(abs(Y-Y'));
    R = radii*ones(1,numberOfAtoms);
    R = triu(R+R',1);
    distance = sqrt(X.^2+Y.^2);
    overlap = distance-R < 0;
    overlapPositions = find(overlap == true);
    for n=1:size(overlapPositions)
        i = mod(overlapPositions(n),numberOfAtoms);
        j = ceil(overlapPositions(n)/numberOfAtoms);
        [velocities(i,:),velocities(j,:)] = impulseFunc(coordinates(i,:),coordinates(j,:),...
            velocities(i,:),velocities(j,:));
    end
    
    %      for i=1:numberOfAtoms-1
    %         for j=i+1:numberOfAtoms
    %             if overlapCheck(coordinates(i,:),coordinates(j,:),radii(i),radii(j))
    %                 %coordinates([i j],:) = coordinates([i j],:) - velocities([i j],:)*dt; %rewind overlapping of cetain atoms
    %                 [velocities(i,:),velocities(j,:)] = impulseFunc(coordinates(i,:),coordinates(j,:),...
    %                     velocities(i,:),velocities(j,:));
    %                 %coordinates([i j],:) = coordinates([i j],:) + velocities([i j],:)*dt; %update
    %             end
    %         end
    %     end
    %
    
    %keep them inside the box
    mirrorVelocities = -([coordinates(:,1)+radii > boxSize(1)...
        coordinates(:,2)+radii > boxSize(2)] +...
        (coordinates-[radii radii] < 0));  %get matrix containing -1 if outside the box
    velocities = (mirrorVelocities+abs(mirrorVelocities+1)).*velocities; %mirror the affected velocity directions
    
    coordinates = coordinates + velocities*dt;
    
    
    if  sum(sum(velocities > (radius/dt))) >= 1
        disp(velocities(find(sqrt(velocities(:,1).^2+velocities(:,2).^2) > (radius/dt) == 1,1),:));
        %error('Out of control')
    end
    
    
    
    if mod(t,1) == 0
        leftSideIndex = find(leftSide==1);
        rightSideIndex = find(leftSide==0);
        viscircles(coordinates(leftSideIndex,:),radii(leftSideIndex),'EdgeColor','b');
        viscircles(coordinates(rightSideIndex,:),radii(rightSideIndex),'EdgeColor','r');
        
        pause(0.02)
    end
    
end