function COE = rvECI2coe(r_vect, v_vect, mu)
% Description:
% Evaluates the classical orbital parameters @ a given epoch starting 
% from the position and velocity vectors.

% Input: 
% r_vect contains (x, y, z) coordinates @ t0
% v_vect contains (vx, vy, vz) @ t0

% Output:
% COE vector contains the Classical Orbital Elements in the following order
% COE = [a, e, incl, Omega, omega, nu]


% First, check the size of position and velocity vectors
if size(r_vect) ~= [3, 1] | size(v_vect) ~= [3, 1]
    fprintf('Error! Dimensions of the vectors provided in rvECI2coe are not correct!')
end

% Compute Auxiliary Parameters
h_vect = cross(r_vect, v_vect);
n_vect = cross([0 0 1], h_vect);    % direction of the line of nodes
n = norm(n_vect); 
r = norm(r_vect);
e_vect = (1/mu) * cross(v_vect, h_vect) - r_vect/r;

% Compute Classical Orbital Parameters
h = norm(h_vect);       % angular momentum
p = h^2 / mu;          % semi-latus rectum
e = norm(e_vect);       % eccentricity 
a = p / (1-e^2);        % semi-major axis

h_vers = h_vect / h;

incl = acos(dot([0 0 1], h_vers));


% Assing values to Omega, omega and nu -> we must follow this procedure
% to solve the ambiguity:

Omega = 0;
omega = 0;
nu = 0;

% Omega > pi IF n_vect(2) < 0
if n_vect(2) >= 0 
    Omega = acos(n_vect(1) / n);
elseif n_vect(2) < 0
    Omega = 2*pi - acos(n_vect(1) / n);
end

% omega > pi IF e_vect(3) < 0
if e_vect(3) >= 0
    omega = acos(dot(n_vect,e_vect) / (n * e));
elseif e_vect(3) < 0
    omega = 2*pi - acos(dot(n_vect,e_vect) / (n * e));
end

% nu > 0 IF dot(r_vect, v_vect) > 0
if dot(r_vect, v_vect) >= 0
    nu = real(acos(dot(e_vect, r_vect) / (e * r)));
elseif dot(r_vect,v_vect) < 0
    nu = 2*pi - real(acos(dot(e_vect, r_vect) / (e * r)));
end


% Value assignment in COE vector
COE = [a, e, incl, Omega, omega, nu];

end


