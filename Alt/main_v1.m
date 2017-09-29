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

numberOfAtoms = 20;             %Number of Atoms
radius = 71;      %Radius of Atoms
itMax = 10000;      %Maximum iterations for spawning atoms
dt = 1;          %Time
velocity = 10;     %Velocity
acceleration = 0;   %Acceleration - will be calculated later via F = m*a


%% jetzt gehts los

radii = ones(numberOfAtoms,1)*radius;

%% Spawning inside box without overlaping box margins
coordinates = [(rand(numberOfAtoms,1)*(boxSize(1)-2*radius))+radius...
    (rand(numberOfAtoms,1)*(boxSize(2)-2*radius))+radius];

%% setting random directions and apply velocity
direction = rand(numberOfAtoms,2)-0.5; %assign random directions
directionNorm = sqrt(direction(:,1).^2 + direction(:,2).^2); %calc normVector
direction = direction./[directionNorm directionNorm]; %normalizing direction
velocities = velocity*direction; %assign velocity

%% let the atoms fly
while 1
    cla
    coordinates = coordinates + velocities*dt;
    
    %keep them inside the box
    mirrorVelocities = -([coordinates(:,1)+radii > boxSize(1)...
        coordinates(:,2)+radii > boxSize(2)] +...
        (coordinates-[radii radii] < 0));  %get matrix containing -1 if outside the box
    velocities = (mirrorVelocities+abs(mirrorVelocities+1)).*velocities; %mirror the affected velocity directions
    
    
    
    viscircles(coordinates,...
        radii);
    
    pause(0.02)
end