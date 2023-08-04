function H = H_matrix_Jacobian(X_sat, Rgs)
% Description: this function uses the Jacobian and Syms functions
% to compute the H matrix for the problem.

% Potrebbe essere sbagliata perchè non differentiable

global log j

syms x y z vx vy vz

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
rho = [x y z]' - Rgs;
rho_hat = rho/norm(rho);

c3 = [0 0 1]';

% Compute Topocentric Declination
Dec = (pi/2 - acos(dot(c3,rho_hat)));           % rad

% Projection on the (x,y) Plane
rho_xy = rho - [0;0;rho(3)];
rho_xy_hat = rho_xy/norm(rho_xy);

% Compute Topocentric Right Ascension
Ra = (atan2(rho_xy_hat(2), rho_xy_hat(1)));     % rad

% Phase Unwrapping on Ra
Ra_val = double(subs(Ra, [x y z vx vy vz], X_sat'));

if Ra_val < 0
    fprintf('Ra is negative at iteration %.0f\n', j)
    Ra_val = Ra_val + 2*pi;
    Ra = Ra + 2*pi;
end


% Assign the Values
Ra_J = jacobian(Ra, [x y z vx vy vz]);
Dec_J = jacobian(Dec, [x y z vx vy vz]);

Ra_J = double(subs(Ra_J, [x y z vx vy vz], X_sat'));
Dec_J = double(subs(Dec_J, [x y z vx vy vz], X_sat'));

H(1, :) = Ra_J;
H(2, :) = Dec_J;

fprintf(log, 'Ra_val = %.6f\n', Ra_val);


end