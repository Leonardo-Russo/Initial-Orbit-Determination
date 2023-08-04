function H = H_matrix(X_sat, Rgs)
% Description: this function uses differential corrections to compute 
% the H matrix for the problem.

global log

% Check Dimensions
if size(Rgs, 2) ~= 1
    Rgs = Rgs';
end

if size(X_sat, 2) ~= 1
    X_sat = X_sat';
end

% Read Input Data
xgs = Rgs(1);
ygs = Rgs(2);
zgs = Rgs(3);

% Initialize H matrix
H = zeros(2, 6);

% Satellite Position wrt Ground Station
rho = X_sat(1:3,1) - Rgs;
rho_hat = rho/norm(rho);

c3 = [0 0 1]';

% Compute Topocentric Declination
Dec = (pi/2 - acos(dot(c3,rho_hat)));           % rad

% Projection on the (x,y) Plane
rho_xy = rho - [0;0;rho(3)];
rho_xy_hat = rho_xy/norm(rho_xy);

% Compute Topocentric Right Ascension
Ra = atan2(rho_xy_hat(2), rho_xy_hat(1));     % rad

% Phase Unwrapping on Ra
if Ra < 0
    Ra = Ra + 2*pi;
end


delta = 1e-10;      % finite difference step size

% Assign the Values
for i = 1 : 6
    dX = zeros(6,1);
    dX(i) = delta;
    
    X_plus = X_sat + dX;
    X_minus = X_sat - dX;

    % Calculate forward difference
    Ra_plus = observation_model(X_plus, Rgs, true);
    Ra_minus = observation_model(X_minus, Rgs, true);
    H(1,i) = (Ra_plus - Ra_minus) / (2*delta);

    % Calculate backward difference
    Dec_plus = observation_model(X_plus, Rgs, false);
    Dec_minus = observation_model(X_minus, Rgs, false);
    H(2,i) = (Dec_plus - Dec_minus) / (2*delta);
end

fprintf(log, 'Ra_val = %.6f\n', Ra);

end

function obs = observation_model(X, Rgs, isRa)
    % Satellite Position wrt Ground Station
    rho = X(1:3,1) - Rgs;
    rho_hat = rho/norm(rho);

    c3 = [0 0 1]';

    % Compute Topocentric Declination
    Dec = (pi/2 - acos(dot(c3,rho_hat)));           % rad

    % Projection on the (x,y) Plane
    rho_xy = rho - [0;0;rho(3)];
    rho_xy_hat = rho_xy/norm(rho_xy);

    % Compute Topocentric Right Ascension
    Ra = atan2(rho_xy_hat(2), rho_xy_hat(1));     % rad

    % Phase Unwrapping on Ra
    if Ra < 0
        Ra = Ra + 2*pi;
    end

    if isRa
        obs = Ra;
    else
        obs = Dec;
    end
end
