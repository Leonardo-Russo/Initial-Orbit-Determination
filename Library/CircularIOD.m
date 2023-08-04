function [r1_vect, v1_vect] = CircularIOD(X, times, Rgs1, Rgs2, mu, rho1_0)
% Description: this function performs IOD with two optical observations 
% (4 angles), by using the hypothesis of circular orbit.
% The first guess value of rho is rho_0.
% The outputs are Position and Velocity at time t1.

% Retrieve Data from Input
RA1 = X(1, 1);
DEC1 = X(1, 2);
RA2 = X(2, 1);
DEC2 = X(2, 2);

global log

t1 = times(1);
t2 = times(2);

% Compute Observations Position Vectors
L1 = [cos(DEC1)*cos(RA1)
      cos(DEC1)*sin(RA1)
      sin(DEC1)];

L2 = [cos(DEC2)*cos(RA2)
      cos(DEC2)*sin(RA2)
      sin(DEC2)];

% % Old Optimization Approach
% % Define Options for the Minimization
% options = optimset('TolFun', 1e-3, 'Disp', 'Iter');
% 
% % Optimization Process
% minFunction = @(rho1) rho_optimal(rho1, L1, L2, Rgs1, Rgs2, t1, t2);
% rho1 = fminsearch(minFunction, rho1_0, options);

% Iterative Computation of rho1
opt_flag = 1;
iteration = 1;
rho_step = 0.1;             % km
rho1_max = 100000;          % km
J_tol = 1e-6;               % rad/s
imag_tol = 1e-6;

rho1 = rho1_0;              % km - Initial Guess
fprintf(log, 'Debug - Circular Assumption rho1 Iterations:\n');

while opt_flag

    fprintf(log, 'rho1 = %.4f\n', rho1);

    if iteration > 1
        rho1 = rho1 + rho_step;
    end

    if rho1 > rho1_max
        error('No solution was found for rho1, please increase the Tolerance.')
    end
    
    r1_vect = rho1*L1 + Rgs1;
    r = norm(r1_vect);

    % Solve the equation for rho2   
    B = 2*dot(L2, Rgs2);
    C = norm(Rgs2)^2 - r^2;

    rho2 = (-B + sqrt(B^2 - 4*C))/2;
    
    if imag(rho2) > imag_tol
        fprintf('rho1 = %.1f + i%.1f\n', real(rho2), imag(rho2))
        error('The solution of rho2 is Complex!')
    end
    
    r2_vect = rho2*L2 + Rgs2;
    
    % Angular Difference between r1 and r2
    dtheta = acos(dot(r1_vect,r2_vect) / (norm(r1_vect)*norm(r2_vect)));  % rad
    
    dt = seconds(t2-t1);
    
    omega_act = abs(dtheta)/dt;
    omega_th = sqrt(mu/r^3);
    
    J = abs(omega_act - omega_th);

    if J < J_tol
        opt_flag = 0;
    else
        rho1 = rho1 + rho_step;
        iteration = iteration + 1;
    end

end


v1 = omega_act*r;
r1_hat = r1_vect / norm(r1_vect);

% We can evaluate the orbital plane because both r1, r2 belong to it!
h_hat = cross(r1_vect, r2_vect) / norm(cross(r1_vect, r2_vect));

% Finally, we can compute the velocity vector
v1_hat = cross(h_hat, r1_hat) / norm(cross(h_hat, r1_hat));
v1_vect = v1 * v1_hat;


end
