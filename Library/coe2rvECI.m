function [r_vect, v_vect] = coe2rvECI(coe, amu)
% Description:
% Transformation from Classical Orbital Elements back
% into Carthesian Coordinates

% Input: 
% COE vector contains the Classical Orbital Elements in the following order
% COE = [a, e, incl, Omega, omega, nu]

% Output:
% r_vect contains (x, y, z) coordinates @ t1
% v_vect contains (vx, vy, vz) @ t1

%% Define r_vect and v_vect in Perifocal Reference Frame

% assume the two main directions to be
% P = [1 0 0]';
% Q = [0 1 0]';       

% Read values from input
a = coe(1);
e = coe(2);
incl = coe(3);
Omega = coe(4);
omega = coe(5);
nu = coe(6);

p = a * (1-e^2);
r = p / (1 + e*cos(nu));
v_r = sqrt(amu/p) * e*sin(nu);
v_t = sqrt(amu/p) * (1 + e*cos(nu));

r_vectP = [r*cos(nu), r*sin(nu), 0]';
v_vectP = [v_r*cos(nu)-v_t*sin(nu), v_t*cos(nu)+v_r*sin(nu), 0]';

%% Definition of Rotation Matrix and Conversion in ECI Frame
% We can split the Rotation Matrix R into three Rotation
% Matrices R1, R2 and R3.

R1 = [[cos(Omega), -sin(Omega), 0] ;...
     [sin(Omega), cos(Omega), 0] ;...
     [0, 0, 1]];

R2 = [[1, 0, 0];...
     [0, cos(incl), -sin(incl)] ;...
     [0, sin(incl), cos(incl)]];

R3 = [[cos(omega), -sin(omega), 0];...
     [sin(omega), cos(omega), 0];...
     [0, 0, 1]];

% Then, we can build the Perifocal2ECI Rotation Matrix as

R = R1 * R2 * R3;

% and then rotate the r and v vectors

r_vect = R * r_vectP;
v_vect = R * v_vectP;

end
