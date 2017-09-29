%% Author Enrico Scherlies

%% Initialization
clear all;
close all;
clc;
clf;
cla;

%% Box Size and Axis Initialization

boxSize =  1e4*[1 1];   %faktor * [width hight] of the Box

%figure settings
axis equal;
axis([0 boxSize(1) 0 boxSize(2)]);
hold on

%% Declaring and Initializing Variables


density = 1 * 29.97/1e5;    %[pm] per atom (argon)
numberOfAtoms =  floor(density^2*boxSize(1)*boxSize(2)); %Calculate Number of Atoms according to Boxsize and Density
disp(['Number of Atoms =' num2str(numberOfAtoms)]);
radius = 71;                %Radius of Atoms in pm
velocity = 416e12;          %Velocity in pm/s
accelerations = zeros(numberOfAtoms,2);     %Acceleration - will be calculated later via F = m*a
radii = ones(numberOfAtoms,1)*radius;       %Radii for each atom

dt = 40e-15; %[s]

%% Lennard-Jones-Potential Values and Constants
sigma = 340.5;          %[pm]
sigma6 = sigma^6;       %used to make computation faster
sigma12 = sigma6^2;     %-||-
epsilon = 1.653e-21;    %10-21 J
mass = 6.634e-26*1e12;  %kg -> pg  (kilograms) (39.948u)


%% PART I: Spawning inside box without overlaping box margins

coordinates = zeros(numberOfAtoms,2);                       %precreate coordinates matrix
coordinates(1,:) = [(rand*(boxSize(1)-2*radius))+radius...  %asign first atom
    (rand*(boxSize(2)-2*radius))+radius];

for i=2:numberOfAtoms %spawn other atoms while checking if newly spawned atom overlaps with any otehr atom

    coordinates(i,:) = [(rand*(boxSize(1)-2*radius))+radius...    %create coordinates:�
        (rand*(boxSize(2)-2*radius))+radius];
    
    %check wether newly created atom overlaps:
    check = false;
    checkIterations = 0;
    while ~check
        checkIterations = checkIterations + 1;
        for j=1:i-1
            if overlapCheck(coordinates(i,:),coordinates(j,:),radii(i),radii(j)) %funkction checking wether two atoms overlap
                coordinates(i,:) = [(rand*(boxSize(1)-2*radius))+radius...
                    (rand*(boxSize(2)-2*radius))+radius];
                break
            end
            if j==i-1
                check = true;
            end
        end
        if checkIterations > 10000  %after a certain numbre of iterations an error message will be displayed
            error('Spawning Atoms without overlap might not possible')
        end
    end
    
end

%% setting random directions and apply velocity
direction = rand(numberOfAtoms,2)-0.5;                       %assign random directions
directionNorm = sqrt(direction(:,1).^2 + direction(:,2).^2); %calc normVector
direction = direction./[directionNorm directionNorm];        %normalizing direction
velocities = velocity*direction;                             %assign velocity

%% Save LeftSide Atom Positions and Draw Plot
leftSide = coordinates(:,1) < boxSize(1)/2; %Saves Positions of Atoms in the Left Side
figure(1)
leftSideIndex = find(leftSide==1);      %Saves the index of the positions either for left and
rightSideIndex = find(leftSide==0);     %right side
viscircles(coordinates(leftSideIndex,:),radii(leftSideIndex),'EdgeColor','g');   %plots left side green
viscircles(coordinates(rightSideIndex,:),radii(rightSideIndex),'EdgeColor','b'); %plots right side blue

disp('Programm paused. Press any key to continue');
pause


%% PART II: Let the atoms fly!
disp('Programm will continue until termination by use (Ctrl + C)');

%Preallocating some Variables:
accelerationArrayX = zeros(numberOfAtoms);
accelerationArrayY = zeros(numberOfAtoms);

velocityDiff = 0;
velocitySumStart = sum(sqrt(velocities(:,1).^2+velocities(:,2).^2)); %== velocity*numberOfAtoms

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
    %disp(sum(sum(velocities)));
    velocitySum = 0;
    velocityDiff = 0;
    
    %calculate overlapping atoms
    X = coordinates(:,1)*ones(1,numberOfAtoms); %get X coordinates
    X = triu(abs(X-X'));                        %get all delta_x_ for every atom 
    Y = coordinates(:,2)*ones(1,numberOfAtoms); %get Y coordiantes
    Y = triu(abs(Y-Y'));                        %calc every delta y
    R = radii*ones(1,numberOfAtoms);            %get Radii
    R = triu(R+R',1);                           %calc sum of radius pairs
    distance = sqrt(X.^2+Y.^2);                 %calc distance between every atom
    overlap = triu(distance-R <= 0,1);          %if distance is less than the radius, overlap becomes true
    overlapPositions = find(overlap == true);   %find overlap positions
    
    %calculate interaction forces of atoms inside 6-sigma intervall
    for n=1:size(overlapPositions)
        i = mod(overlapPositions(n),numberOfAtoms);  %get atom positions i and
        j = ceil(overlapPositions(n)/numberOfAtoms); %j
        
        [accAtom1,accAtom2] =  accLJP(coordinates(i,:),coordinates(j,:),sigma6,sigma12,epsilon,mass); %use accLJP Function
        accelerationArrayX(i,j) = accAtom1(1); %save in X- and Y-Arrays
        accelerationArrayX(j,i) = accAtom2(1);
        accelerationArrayY(i,j) = accAtom1(2);
        accelerationArrayY(j,i) = accAtom2(2);
    end

    accelerations = [sum(accelerationArrayX)' sum(accelerationArrayY)']; %sum all vectors to get new acceleration
    
    %update new coordinates and velocities
    velocities = velocities + accelerations*dt;
    coordinates = coordinates + velocities*dt + 0.5*accelerations*(dt.^2);
    
    
    %Energy conservation
%      velocitySum = sum(sqrt(velocities(:,1).^2+velocities(:,2).^2));
%     velocityDiff = (velocitySumStart - velocitySum)/numberOfAtoms;
%     dn(velocityDiff);
%     velocities = velocities + ones(numberOfAtoms,2)*velocityDiff;
    
    
    %Keep atoms inside the box
    mirrorVelocities = -([coordinates(:,1)+radii > boxSize(1)...
        coordinates(:,2)+radii > boxSize(2)] +...
        (coordinates-[radii radii] < 0));  %get matrix containing -1 if outside the box
    velocities = (mirrorVelocities+abs(mirrorVelocities+1)).*velocities; %mirror the affected velocity directions
    
    
    
    %dt = radius/max(sqrt(velocities(:,1).^2+velocities(:,2).^2));
    
    
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
        %pause(0.02)
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