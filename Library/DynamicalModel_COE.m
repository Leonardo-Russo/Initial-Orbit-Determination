function dx = DynamicalModel_COE(t, X)
% Description: This will be the function to integrate to obtain the
% propagation.
% Input: t = tspan(i), X = COE = [a, e, incl, Omega, omega, nu]

% Recalling Global Variables
global mu

% Reading State Values
N = length(X);

a = X(1);
e = X(2);
incl = X(3);
Omega = X(4);
omega = X(5);
nu = X(6);

fr = 0;
ft = 0;
fh = 0;

coe = [a, e, incl, Omega, omega, nu];

% Initialize State Derivative
dx = zeros(N, 1);

% Introduce Local Variables
theta = omega + nu;
p = a*(1-e^2);

[r_vect, v_vect] = coe2rvECI(coe, mu);
r = norm(r_vect);
h_vect = cross(r_vect, v_vect);
h = norm(h_vect);


% Assign Values to State Derivative
dx(1) = 2*a^2/h * (e*fr*sin(nu) + ft*(e*cos(nu) + 1));
dx(2) = sqrt(p/mu)*fr*sin(nu) + sqrt(p/mu)*ft*(e+e*cos(nu)^2+2*cos(nu))/(1+e*cos(nu));
dx(3) = r*fh/h * cos(theta);
dx(4) = r*fh/h *sin(theta)/sin(incl);
dx(5) = -r*fh*sin(theta)*cos(incl)/(h*sin(incl)) - fr/e*cos(nu)*sqrt(p/mu) + ft*sqrt(p/mu)*sin(nu)*(e*cos(nu)+2)/(e*(1+e*cos(nu)));
dx(6) = sqrt(mu/p^3)*(1+e*cos(nu))^2 + fr/e*cos(nu)*sqrt(p/mu) - ft*sqrt(p/mu)*sin(nu)*(e*cos(nu)+2)/e*(1+e*cos(nu));


end