%% A2 EOSC 453 
% Max, Joel and Keegan

%% Define Input Parameters
clear; close all; clear global
% Make variables that will be used in the ODE global
global k L ao a pcZaverage sigmaB S0 tau asky gamma A phi1 phi2 ...
        phi3 phi4 phi5 phi6 ai Tlow Thigh count;

% Load RNG Seed
load('seed1.mat');

% Reset RNG
rng(s);

% Thermal Exchange Coefficients (W/mK)
k = zeros(6);
k(1,2) = 1*10^7;
k(2,3) = 1*10^7;
k(3,4) = 1*10^7;
k(4,5) = 5*10^7;
k(5,6) = 1*10^7;

% Boundary Length (m)
L = zeros(6);
L(1,2) = 2.0015*10^7;
L(2,3) = 3.4667*10^7;
L(3,4) = 4.0030*10^7;
L(4,5) = 3.4667*10^7;
L(5,6) = 2.0015*10^7;

% Area Fractions of Zones
areafrac = [0.0670 0.1830 0.2500 0.2500 0.1830 0.0670];

% Geometric Factor
gamma = [0.1076 0.2277 0.3045 0.3045 0.2277 0.1076];

% Constants
sigmaB = 5.6696*10^-8; % W m?2 K?4 Stefan-Boltzmann constant
S0 = 1368; % W m?2 Solar constant
RE = 6371*10^3; %m Radius of Earth, 
AE = 4*pi*RE^2; % Surface area of Earth
eta = 1; % Total emissivity of Earth
tau = 0.63; % Atmospheric transmissivity 
asky = 0.2; % Atmospheric albedo
al = 0.4; % Albedo of land surface
as = 0.1; % Albedo of ocean surface
ai = 0.6; % Albedo of ice
pl = 2500; % kg m?3 Density of land surface
ps = 1028; % kg m?3 Density of ocean water
pi = 900; % kg m?3 Density of ice
Zl = 1.0; % m Thermal scale depth for land
Zs = 70.0; % m Thermal scale depth for ocean
Zi = 1.0; % m Thermal scale depth for ice
cl = 790; % J kg K?1 Specific heat capacity for land
cs = 4187; % J kg K?1 Specific heat capacity for water
ci = 2060; % J kg K?1 Specific heat capacity for ice
fl = [0.05 0.25 0.3 0.45 0.55 0.35]; % Land Fraction
fs = [0.54 0.74 0.7 0.55 0.44 0.45]; % Sea Fraction
fi = [0.41 0.01 0.0 0.00 0.01 0.20]; % Ice Fraction
sec_per_year = 365*24*3600; % Seconds per year
Tlow = 260;
Thigh = 273;
count = 1;

% Area of Zones
A = AE*areafrac;

% Initial Temperature Conditions
To = [290 290 290 290 290 290];

% Time Parameters
starttime_Ma = 0;
starttime_s = starttime_Ma*1000000*365*24*3600;

endtime_Ma = 0.01;
endtime_s = endtime_Ma*1000000*365*24*3600;

% Surface Average Values

pcZaverage(1)=fl(1)*pl*cl*Zl+fs(1)*ps*cs*Zs+fi(1)*pi*ci*Zi;
pcZaverage(2)=fl(2)*pl*cl*Zl+fs(2)*ps*cs*Zs+fi(2)*pi*ci*Zi;
pcZaverage(3)=fl(3)*pl*cl*Zl+fs(3)*ps*cs*Zs+fi(3)*pi*ci*Zi;
pcZaverage(4)=fl(4)*pl*cl*Zl+fs(4)*ps*cs*Zs+fi(4)*pi*ci*Zi;
pcZaverage(5)=fl(5)*pl*cl*Zl+fs(5)*ps*cs*Zs+fi(5)*pi*ci*Zi;
pcZaverage(6)=fl(6)*pl*cl*Zl+fs(6)*ps*cs*Zs+fi(6)*pi*ci*Zi;

% Avg Albedo
ao(1)= fl(1)*al+fs(1)*as+fi(1)*ai;
ao(2)= fl(2)*al+fs(2)*as+fi(2)*ai;
ao(3)= fl(3)*al+fs(3)*as+fi(3)*ai;
ao(4)= fl(4)*al+fs(4)*as+fi(4)*ai;
ao(5)= fl(5)*al+fs(5)*as+fi(5)*ai;
ao(6)= fl(6)*al+fs(6)*as+fi(6)*ai;
a = zeros(size(ao));
%% Calc ODE for Steady State
[time, Temps] = ode15s(@tempodes,[starttime_s endtime_s], To);

figure(1)
plot(time/sec_per_year,Temps)
legend('1','2','3','4','5','6')
title('Steady State (ODE15s)')

t1 = Temps(end,:);

%% Random Events Effecting the Incoming Solar Radiation Coefficient (Phi)
% We want it to have random magnitudes, at random times, in random zones.
% This section will initially create a phi that is all ones, so that it 
% will not effect steady state. 

% Number of Volcanic Events
num_events = 100;

% Create Incoming Solar Radiation Coefficient for each zone
phi1(:,1) = starttime_s:sec_per_year:endtime_s; %s
phi1(:,2) = ones(size(phi1(:,1)));
phi2=phi1;
phi3=phi1;
phi4=phi1;
phi5=phi1;
phi6=phi1;

% Find a Zone
zone_event = randi([1 100],1,num_events);

% Find percent of sphere each zone takes up
areaperc = round(areafrac*100);

% Distribute zones according to percent of sphere
for g = 1:length(areaperc);
 
   gg = zone_event <= sum(areaperc(1:g)) & zone_event > sum(areaperc(1:g-1));
   zone_event(gg) = g;
   
end
% Get a Magnitude
mag_event = 0.5*rand(1,num_events);

%Get a point in time
time_event = randi([starttime_s endtime_s]/sec_per_year,1,num_events); %s

% Set a linear rebound time scale for phi after an event. This value will be
% multiplied with the magnitude to determine how many years it takes to
% return to normal
rebound = 10;

% Loop over the number of volcanic events
for i = 1:num_events;

    %Find where it is in phi's time
    [dummy, closest_time] = min(abs(phi1(:,1) - time_event(i)*sec_per_year));

    %Define event
    event = linspace((1 - mag_event(i)), 1, rebound*ceil(mag_event(i)));
    
    % Index the rows in phi that will be replaced with the event
    ii = closest_time:closest_time+length(event)-1;

    % Replace values in phi with the event and distribute to zones
    switch zone_event(i)
        
        case 1    
            phi1(ii,2) = event;
        case 2    
            phi2(ii,2) = event;
        case 3    
            phi3(ii,2) = event;
        case 4    
            phi4(ii,2) = event;
        case 5    
            phi5(ii,2) = event;
        case 6    
            phi6(ii,2) = event;
    end
    
end

% %% Calc ODE for Phi Forcing
% [time_phi, Temps_phi] = ode45(@tempodes_phi,[starttime_s endtime_s], To);
% 
% figure(2)
% plot(time_phi/sec_per_year,Temps_phi)
% legend('1','2','3','4','5','6')
% 
% figure
% scatter(time_event,zone_event)

%% Calc ODE for Phi Forcing w/ Ice - Albedo Feedback
[time_phi_AIF, Temps_phi_AIF] = ode45(@tempodes_phi_feedback,[starttime_s endtime_s], To);

figure(3)
plot(time_phi_AIF/sec_per_year,Temps_phi_AIF)
legend('1','2','3','4','5','6')





