function [r_vect, v_vect] = LaplaceIOD(X, times, Rgs, dRgs, ddRgs, t, mu)
% Description: this function implements the Laplace Method for Initial
% Orbit Determination.
% X, times = Observations State and associated Observation Times
% Rgs, dRgs, ddRgs = Position, Velocity, Acceleration of GS in ECI
% t = time at which to perform the IOD method
% P = constant parameters
% r_vect, v_vect = Position, Velocity of the S/C @ t.


% Retrieve Data from Input
RA1 = X(1, 1); 
DEC1 = X(1, 2);
RA2 = X(2, 1);
DEC2 = X(2, 2);
RA3 = X(3, 1);
DEC3 = X(3, 2);

t1 = times(1);
t2 = times(2);
t3 = times(3);


% Compute Observations Position Vectors
L1 = [cos(DEC1)*cos(RA1)
      cos(DEC1)*sin(RA1)
      sin(DEC1)];

L2 = [cos(DEC2)*cos(RA2)
      cos(DEC2)*sin(RA2)
      sin(DEC2)];

L3 = [cos(DEC3)*cos(RA3)
      cos(DEC3)*sin(RA3)
      sin(DEC3)];


% Rewrite times wrt a reference time t0
t0 = datetime([2022, 1, 1, 0, 0, 0]);
ts = seconds(t-t0);
t1s = seconds(t1-t0);
t2s = seconds(t2-t0);
t3s = seconds(t3-t0);

% Compute Satellite Position Versor L via Interpolation
L = (ts-t2s)*(ts-t3s)/((t1s-t2s)*(t1s-t3s))*L1 ...
    + (ts-t1s)*(ts-t3s)/((t2s-t1s)*(t2s-t3s))*L2 ...
    + (ts-t1s)*(ts-t2s)/((t3s-t1s)*(t3s-t2s))*L3;

% Angular Velocity of L
dL = (2*ts-t2s-t3s)/((t1s-t2s)*(t1s-t3s))*L1 ...
    + (2*ts-t1s-t3s)/((t2s-t1s)*(t2s-t3s))*L2 ...
    + (2*ts-t1s-t2s)/((t3s-t1s)*(t3s-t2s))*L3;

% Angular Acceleration of L
ddL = 2/((t1s-t2s)*(t1s-t3s))*L1 ...
    + 2/((t2s-t1s)*(t2s-t3s))*L2 ...
    + 2/((t3s-t1s)*(t3s-t2s))*L3;


% Compute Local Variables
D = 2*det([L, dL, ddL]);
D1 = det([L, dL, ddRgs]);
D2 = det([L, dL, Rgs]);

% Solve the Equation for rho and r
syms rho r

eq1 = rho + 2*D1/D + 2*mu/r^3*D2/D;
eq2 = r^2 - (rho^2 + 2*rho*dot(L, Rgs) + dot(Rgs, Rgs));
sol = solve([eq1, eq2], [rho, r]);

n_sol = 0;
rhos = [];
rs = [];

for i = 1 : length(sol.rho)

    rhos = [rhos; double(sol.rho(i))];
    rs = [rs; double(sol.r(i))];

end

% Only select Real Solutions
for i = 1 : length(rhos)
    if real(rhos(i)) > 0 && imag(rhos(i)) == 0 && real(rs(i)) > 0 && imag(rs(i)) == 0
        rho = real(rhos(i));
        r = real(rs(i));
        n_sol = n_sol + 1;
    end
end

if n_sol > 1
    warning("There's more than one acceptable solution!")
end

if n_sol < 1
    error("There's no acceptable solution!")
end

r_vect = Rgs + rho*L;

D3 = det([L, ddRgs, ddL]);
D4 = det([L, Rgs, ddL]);

drho = -D3/D - mu/r^3 * D4/D;

v_vect = drho*L + rho*dL + dRgs;


end