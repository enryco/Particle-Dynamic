%% Initialization
clear all;
close all;
clc;
clf;
cla;

%% Box Size and Axis Initialization

boxSize = 8e4*[1 1];   %[width hight] of the Box

axis equal;
axis([0 boxSize(1) 0 boxSize(2)]);
hold on

%% Declaring Variables

density = 29.97/1e5; %pm per atom (argon)
numberOfAtoms =  floor(density^2*boxSize(1)*boxSize(2)); %Number of Atoms
disp(['Number of Atoms =' num2str(numberOfAtoms)]);
radius = 71;        %Radius of Atoms in pm
itMax = 10000;      %Maximum iterations for spawning atoms
%dt = 0.1;             %Time
velocity = 416*1e12;      %Velocity in pm/s
accelerations = zeros(numberOfAtoms,2);   %Acceleration - will be calculated later via F = m*a

dt = radius/velocity/100;

%LJP constants
sigma = 340.5; %pm
sigma6 = sigma^6;
sigma12 = sigma6^2;
epsilon = 1.653e-21; % 10-21 J
mass = 6.634e-26*1e12; %kg -> pg  (kilograms) (39.948u)

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

leftSide = coordinates(:,1) < boxSize(1)/2;
%circleColor = [leftSide zeros(numberOfAtoms,1) not(leftSide)];


%% setting random directions and apply velocity
direction = rand(numberOfAtoms,2)-0.5; %assign random directions
directionNorm = sqrt(direction(:,1).^2 + direction(:,2).^2); %calc normVector
direction = direction./[directionNorm directionNorm]; %normalizing direction
velocities = velocity*direction; %assign velocity


%% test
% coordinates = [[250 500];[750 500]];
% velocities = [[velocity velocity];[-velocity velocity]];

%% draw
figure(1)
leftSideIndex = find(leftSide==1);
rightSideIndex = find(leftSide==0);
viscircles(coordinates(leftSideIndex,:),radii(leftSideIndex),'EdgeColor','b');
viscircles(coordinates(rightSideIndex,:),radii(rightSideIndex),'EdgeColor','g');

pause
%clf
%% declairing some usefull v
accelerationArrayX = zeros(numberOfAtoms);
accelerationArrayY = zeros(numberOfAtoms);

velocityDiff = 0;
velocitySumStart = sum(sqrt(velocities(:,1).^2+velocities(:,2).^2)); %== velocity*numberOfAtoms


%% let the atoms fly
iteration = 0;
iterations = [];
concentrationBlue = [];
concentrationGreen = [];
hold on
%cla
while 1
    iteration = iteration+1;
    %cla
    
    accelerationArrayX = zeros(numberOfAtoms);
    accelerationArrayY = zeros(numberOfAtoms);
    %velocitySum = 0;
    %velocityDiff = 0;
    
    
    %calculate new acceeleratsoins
    for i=1:numberOfAtoms-1
        for j=(i+1):numberOfAtoms
            [accAtom1,accAtom2] =  accLJP(coordinates(i,:),coordinates(j,:),sigma6,sigma12,epsilon,mass);
            accelerationArrayX(i,j) = accAtom1(1);
            accelerationArrayX(j,i) = accAtom2(1);
            accelerationArrayY(i,j) = accAtom1(2);
            accelerationArrayY(j,i) = accAtom2(2);
        end
    end
    
    accelerations = [sum(accelerationArrayX)' sum(accelerationArrayY)'];
    
    
    
    %%update new coordinates and velocities
    velocities = velocities + accelerations*dt;
    coordinates = coordinates + velocities*dt + 0.5*accelerations*(dt.^2);
    
    
    velocitySum = sum(sqrt(velocities(:,1).^2+velocities(:,2).^2));
    
    %Energy conservation
    velocityDiff = (velocitySumStart - velocitySum)/numberOfAtoms;
    %velocities = velocities + ones(numberOfAtoms,2)*velocityDiff;
    
    
    %keep them inside the box
    mirrorVelocities = -([coordinates(:,1)+radii > boxSize(1)...
        coordinates(:,2)+radii > boxSize(2)] +...
        (coordinates-[radii radii] < 0));  %get matrix containing -1 if outside the box
    velocities = (mirrorVelocities+abs(mirrorVelocities+1)).*velocities; %mirror the affected velocity directions
    
    
    
    dt = radius/max(sqrt(velocities(:,1).^2+velocities(:,2).^2));
    
    
    %     %     velocitySumStart/sum(velocities)
    %     if (velocitySumStart/sum(sqrt(velocities(:,1).^2+velocities(:,2).^2)) == not(1)) ||...
    %             sum(sum(velocities > (radius/dt))) >= 1
    %         disp(velocities(find(sqrt(velocities(:,1).^2+velocities(:,2).^2) > (radius/dt) == 1,1),:));
    %         disp(velocitySumStart/sum(sqrt(velocities(:,1).^2+velocities(:,2).^2)));
    %         %error('Out of control')
    %     end
    
    
    %draw concentratoin graph
    if mod(iteration,10) == 0
        figure(2)
        cla
        hold on
        iterations = [iterations iteration];
        leftSideCurrent = coordinates(:,1) < boxSize(1)/2;
        concentrationBlue = [concentrationBlue...
            sum((leftSide + leftSideCurrent) == 2)];
        concentrationGreen = [concentrationGreen...
            sum((not(leftSide) + leftSideCurrent) == 2)];
        plot(iterations,concentrationBlue,'.-b');
        plot(iterations,concentrationGreen,'.-g');
        pause(0.02)
    end
    
    
    %%Draw atoms
    if mod(iteration,1) == 0
        figure(1)
        clf
        axis equal;
        axis([0 boxSize(1) 0 boxSize(2)]);
        leftSideIndex = find(leftSide==1);
        rightSideIndex = find(leftSide==0);
        viscircles(coordinates(leftSideIndex,:),radii(leftSideIndex),'EdgeColor','g');
        viscircles(coordinates(rightSideIndex,:),radii(rightSideIndex),'EdgeColor','b');
        %plot(coordinates(:,1),coordinates(:,2),'.');       
        %plotv(velocities');
        %plotv(accelerations');
        %disp(sum(sum(accelerations)))
        pause(0.02)
    end
    
end