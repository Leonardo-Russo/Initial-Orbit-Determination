function [Rgs, dRgs, ddRgs] = groundStation_Old(Long, LatGD, Alt, t)
% Description: this function obtains the GS position vector and its
% derivatives at the specified time from the input.
% Input: [Geographical Longitude, Geodetical Latitude, Elevation, Time]

global thetaG

D_sid = 86164;      % s
D_sol = 86400;      % s

% Define Earth Parameters
aE = 6378.145;              % km
bE = 6356.785;              % km
eE = 0.081771980616944;     
omegaE = 2*pi/D_sid;        % rad/s

xgs = (aE/sqrt(1-eE^2*sin(LatGD)^2) + Alt)*cos(LatGD);
ygs = (aE*(1-eE^2)/sqrt(1-eE^2*sin(LatGD)^2) + Alt)*sin(LatGD);

% Define Reference Greenwich Angle from 01/01/2022 at 00:00
% https://aa.usno.navy.mil/data/siderealtime
thetaG0 = deg2rad(15*(6 + 42/60 + 31.2087/3600));
t0 = datetime(2022, 1, 1, 0, 0, 0);

dt = seconds(t-t0);
thetaG = thetaG0 + omegaE*dt;
theta = thetaG + Long;

Rgs = [xgs*cos(theta)
     xgs*sin(theta)
     ygs];

dRgs = [-xgs*omegaE*sin(theta)
       xgs*omegaE*cos(theta)
            0];

ddRgs = [-xgs*omegaE^2*cos(theta)
       -xgs*omegaE^2*sin(theta)
            0];


end