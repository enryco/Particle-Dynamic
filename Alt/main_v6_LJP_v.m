%% Initialization
clear all;
close all;
clc;
clf;
cla;

%% Box Size and Axis Initialization

boxSize = 1e3*[1 1];   %[width hight] of the Box

axis equal;
axis([0 boxSize(1) 0 boxSize(2)]);
hold on

%% Declaring Variables

numberOfAtoms = 2; %Number of Atoms
radius = 71;        %Radius of Atoms in pm
itMax = 10000;      %Maximum iterations for spawning atoms
dt = 0.1;             %Time
velocity = 416*1e12;      %Velocity in pm/s
accelerations = zeros(numberOfAtoms,2);   %Acceleration - will be calculated later via F = m*a

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



%% setting random directions and apply velocity
direction = rand(numberOfAtoms,2)-0.5; %assign random directions
directionNorm = sqrt(direction(:,1).^2 + direction(:,2).^2); %calc normVector
direction = direction./[directionNorm directionNorm]; %normalizing direction
velocities = velocity*direction; %assign velocity


%% test 
coordinates = [[201 500];[750 500]];
velocities = [[velocity velocity];[-velocity velocity]];

%% draw
viscircles(coordinates,...
    radii);
pause

%% declairing some usefull v
accelerationArrayX = zeros(numberOfAtoms);
accelerationArrayY = zeros(numberOfAtoms);

velocityDiff = 0;
velocitySumStart = sum(sqrt(velocities(:,1).^2+velocities(:,2).^2)) %== velocity*numberOfAtoms


%% let the atoms fly
t = 0;

while 1
    t = t+1;
    %cla
    
    dt = (radius-1)/max(sqrt(velocities(:,1).^2+velocities(:,2).^2));
    dt = dt/100;
    accelerationArrayX = zeros(numberOfAtoms);
    accelerationArrayY = zeros(numberOfAtoms);
    %velocitySum = 0;
    %velocityDiff = 0;
    
    %keep them inside the box
    mirrorVelocities = -([coordinates(:,1)+radii > boxSize(1)...
        coordinates(:,2)+radii > boxSize(2)] +...
        (coordinates-[radii radii] < 0));  %get matrix containing -1 if outside the box
    velocities = (mirrorVelocities+abs(mirrorVelocities+1)).*velocities; %mirror the affected velocity directions
    
    
    %calculate new acceeleratsoins
    for i=1:numberOfAtoms-1
        for j=(i+1):numberOfAtoms
            [accAtom1,accAtom2] =  accLJP(coordinates(i,:),coordinates(j,:));
            accelerationArrayX(i,j) = accAtom1(1);
            accelerationArrayX(j,i) = accAtom2(1);
            accelerationArrayY(i,j) = accAtom1(2);
            accelerationArrayY(j,i) = accAtom2(2);
        end
    end
    
    accelerations = [sum(accelerationArrayX)' sum(accelerationArrayY)'];
    
    
    
    %%update new coordinates and velocities
    coordinates = coordinates + velocities*dt + 0.5*accelerations*(dt.^2);
    velocities = velocities + accelerations*dt;
    
    
    velocitySum = sum(sqrt(velocities(:,1).^2+velocities(:,2).^2));
   
    %Energy conservation
    velocityDiff = (velocitySumStart - velocitySum)/numberOfAtoms;
    %velocities = velocities + ones(numberOfAtoms,2)*velocityDiff;
    
    
    %     velocitySumStart/sum(velocities)
    if (velocitySumStart/sum(sqrt(velocities(:,1).^2+velocities(:,2).^2)) == not(1)) ||...
            sum(sum(velocities > (radius/dt))) >= 1
        disp(velocities(find(sqrt(velocities(:,1).^2+velocities(:,2).^2) > (radius/dt) == 1,1),:));
        disp(velocitySumStart/sum(sqrt(velocities(:,1).^2+velocities(:,2).^2)));
        %error('Out of control')
    end
    
    
    %%Draw
    if mod(t,10) == 0
        %clf
        %viscircles(coordinates,radii);
        plot(coordinates(:,1),coordinates(:,2),'.');
        %plotv(velocities');
        %plotv(accelerations');
        %disp(sum(sum(accelerations)))
        %toc
        pause(0.02)
        %tic
    end
end