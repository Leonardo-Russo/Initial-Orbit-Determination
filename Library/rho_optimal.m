function J = rho_optimal(rho1, L1, L2, Rgs1, Rgs2, t1, t2)
    
global mu

r1_vect = rho1*L1 + Rgs1;
r = norm(r1_vect);  % circular orbit -> r = const.

% Solve the equation for rho2
syms rho2

eq1 = rho2^2 + dot(Rgs2, Rgs2) + 2*rho2*dot(L2, Rgs2) - r^2;
sol = solve(eq1, rho2);

n_sol = 0;
rho2s = [];

for i = 1 : length(sol)
    rho2s = [rho2s; double(sol(i))];
end

% Only select Real Solutions
for i = 1 : length(rho2s)
    if real(rho2s(i)) > 0 && imag(rho2s(i)) == 0
        rho2s_reals = real(rho2s(i));
        n_sol = n_sol + 1;
    end
end

if n_sol > 1
    warning("There's more than one acceptable solution!")
end

if n_sol < 1
    error("There's no acceptable solution!")
end

% DEVI PRENDERE QUELLO CHE DA IL MINIMO
rho2 = max(rho2s_reals);

r2_vect = rho2*L2 + Rgs2;

% Angular Difference between r1 and r2
dtheta = acos(dot(r1_vect,r2_vect) / (norm(r1_vect)*norm(r2_vect)));  % rad

dt = seconds(t2-t1);

omega_act = abs(dtheta)/dt;
omega_th = sqrt(mu/r^3);

J = abs(omega_act - omega_th);
    
end