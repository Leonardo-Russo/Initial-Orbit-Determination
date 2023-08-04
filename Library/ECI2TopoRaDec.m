function [Ra, Dec] = ECI2TopoRaDec(obs, sat)
% Description: the function takes the ECI positions of the observatory and
% the satellite and yields topocentric Ra and Dec in degrees, as seen from
% the observer.

% Check Dimensions
if size(obs, 2) ~= 1
    obs = obs';
end

if size(sat, 2) ~= 1
    sat = sat';
end

% Satellite Position wrt Observer
obs2sat = sat - obs;
obs2sat_hat = obs2sat/norm(obs2sat);

z = [0 0 1]';

% Compute Topocentric Declination
Dec = (pi/2 - acos(dot(z,obs2sat_hat)))*180/pi;

% Projection on the (x,y) Plane
obs2sat_xy = obs2sat - [0;0;obs2sat(3)];
obs2sat_xy_hat = obs2sat_xy/norm(obs2sat_xy);

% Compute Topocentric Right Ascension
Ra = (atan2(obs2sat_xy_hat(2),obs2sat_xy_hat(1)))*180/pi;

% Check that Ra > 0
if Ra < 0
    Ra = Ra + 360;
end


end
