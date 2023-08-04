function [Rgs, dRgs, ddRgs] = lla2RdRddR(Long, LatGD, Alt, epoch, P)
% Description: this function obtains the GS position vector and its
% derivatives at the specified time from the input.
% Input: [Geographical Longitude, Geodetical Latitude, Elevation, Time, Constants]

global thetaG

D_sid = 86164;      % s

% Define Earth Parameters
aE = 6378.1366;                 % km
bE = 6356.7519;                 % km
eE = sqrt(1 - (bE^2/aE^2)); 
omegaE = 2*pi/D_sid;            % rad/s

% Retrieve Data from Input
DAT = P(1);
ut1_utc = P(2);
X_EOP = P(3);
Y_EOP = P(4);
dX_EOP = P(5);
dY_EOP = P(6);

Rgs = lla2eci([rad2deg(LatGD), rad2deg(Long), Alt*1e3], datevec(epoch), 'IAU-2000/2006', DAT, ut1_utc, [X_EOP, Y_EOP], 'dcip', [dX_EOP, dY_EOP])' * 1e-3;   % km

% Retrieve Data from Rgs
ygs = Rgs(3);
theta = atan2(Rgs(2), Rgs(1));
thetaG = theta - Long;
xgs = Rgs(1)/cos(theta);

% % Old Procedure
% xgs_old = (aE/sqrt(1-eE^2*sin(LatGD)^2) + Alt)*cos(LatGD);
% ygs_old = (aE*(1-eE^2)/sqrt(1-eE^2*sin(LatGD)^2) + Alt)*sin(LatGD);


dRgs = [-xgs*omegaE*sin(theta)
       xgs*omegaE*cos(theta)
            0];


ddRgs = [-xgs*omegaE^2*cos(theta)
       -xgs*omegaE^2*sin(theta)
            0];


end