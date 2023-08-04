function dx = DynamicalModel_Kalman(t, X)
% Description: the function reports the dynamical model for the universal
% gravitational law without perturbances.
% Input: t = tspan(i), X = X(i) = [x, y, z, u, v, z]

N = length(X);

global mu

x = X(1);
y = X(2);
z = X(3);
u = X(4);
v = X(5);
w = X(6);

r = sqrt(x^2 + y^2 + z^2);

phi = X(7:end);
Phi = reshape(phi, 6, 6);
A = A_matrix(X(1:3));
dPhi = A*Phi;

dx = zeros(N, 1);

dx(1) = u;
dx(2) = v;
dx(3) = w;
dx(4) = -mu/r^3 * x;
dx(5) = -mu/r^3 * y;
dx(6) = -mu/r^3 * z;
dx(7:end) = reshape(dPhi, 36, 1);


end