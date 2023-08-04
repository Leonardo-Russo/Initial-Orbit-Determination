%% Initial Orbit Determination - Leonardo Russo 2015563

close all
clear all
clc

addpath('Library/')
addpath('Data/')
addpath('SGP4/')

optional_plots = 1;
optimal_IOD = 0;        % only recommended for Sentinel 3A 
saveplots = 0;

global mu log

log = fopen('Data/log.txt', 'w+');

%% Introduce Known Quantities

choice = input('\n    Please, choose the Dataset:\n    1. Sentinel 3A\n    2. Sentinel 3B\n    3. TDM Beidou\n');

% Read Measurement Error from Header of Fits Images
switch choice

    case 1
        load('Sentinel3A.mat')
        TLE_text = 'Sentinel3A_TLE.txt';
        rho1_0 = 100;           % km
        sigma_meas = deg2rad(8/3600);
        sigma_r = 10;           % km
        sigma_v = 0.01;         % km/s
        % Define the Target Set -> inclination
        target_incl = deg2rad(98.562);      % rad

    case 2
        load('Sentinel3B.mat')
        TLE_text = 'Sentinel3B_TLE.txt';
        rho1_0 = 100;       % km
        sigma_meas = deg2rad(8/3600);
        sigma_r = 10;       % km
        sigma_v = 0.01;     % km/s
        % Define the Target Set -> inclination
        target_incl = deg2rad(98.512);      % rad

    case 3
        load('Beidou.mat')
        TLE_text = 'Beidou_TLE.txt';
        rho1_0 = 30000;       % km
        sigma_meas = deg2rad(2/3600);
        sigma_r = 100;      % km
        sigma_v = 0.1;      % km/s

    otherwise
        error('Invalid User Choice Choice')

end

clc

mu = 398600.4415;           % km^3/s^2
rE = 6378.1363;             % km
D_sid = 86164;              % s
D_sol = 86400;              % s

t1 = tspan_obs(1);

% Geocentric Position of Collepardo Observatory
long_site = deg2rad(13 + 22/60 + 9.84/3600);    % rad - longitude
lat_site = deg2rad(41 + 45/60 + 51.4794/3600);  % rad - latitude
alt_site = 576 * 1e-3;                          % km - altitude
latgd_site = deg2rad(geoc2geod(rad2deg(lat_site), rE*1e3)); % rad - geodetic latitude

% Load TLE Data
[tle1, tle2] = get_TLE(TLE_text);

% Find EOP at Observation Time
obs_start_vect = datevec(tspan_obs(1));
EOPf = findEOPf(obs_start_vect);
MJD = EOPf(4);
X_EOP = arcsec2rad(EOPf(5));
Y_EOP = arcsec2rad(EOPf(6));
ut1_utc = EOPf(7);
lod = EOPf(8);
dpsi = arcsec2rad(EOPf(9));
deps = arcsec2rad(EOPf(10));
dX_EOP = arcsec2rad(EOPf(11));
dY_EOP = arcsec2rad(EOPf(12));
DAT = EOPf(13);

% Defining the options for the Integration
Tol0 = 1e-9;
Tol1 = 1e-11;
MaxStep = 1;
options = odeset('RelTol', Tol0, 'AbsTol',Tol1, 'MaxStep', MaxStep);


%% Perform the SGP4 Propagation for Comparison Purposes

[~, ~, ~, satrec] = twoline2rv(tle1, tle2, 'c', [], [], 84);

% Compute TLE epoch and Check Validity
TLE_epoch = datetime(satrec.jdsatepoch + satrec.jdsatepochf,'ConvertFrom','juliandate');

XECI_SGP4 = zeros(N, 6);
XECEF_SGP4 = zeros(N, 6);

for k = 1 : N   % propagation in time
        
    obs_epoch = tspan_obs(k);
    obs_epoch_vect = datevec(obs_epoch);
    [obs_y, obs_mo, obs_d, obs_h, obs_mi, obs_s] = split_epochvec(obs_epoch_vect);

    % Observatory Coordinates in ECI
    obs_ECI = lla2eci([rad2deg(latgd_site), rad2deg(long_site), alt_site*1e3], obs_epoch_vect, 'IAU-2000/2006', DAT, ut1_utc, [X_EOP, Y_EOP], 'dcip', [dX_EOP, dY_EOP])' * 1e-3;   % km

    % Compute tsince as time elapsed from TLE to Observation in minutes
    tsince = minutes(obs_epoch - TLE_epoch);

    % Propagation of position to Observation Time
    [~, r_vectTEME, v_vectTEME] = sgp4(satrec, tsince);
    
    % Evaluate ttt from UTC
    [~, ~, jdut1, jdut1f, ~, ~, ~, ttt, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~]...
        = convtime(obs_y, obs_mo, obs_d, obs_h, obs_mi, obs_s, 0, ut1_utc, DAT);     % timezone has been set to 0
    
    jdut1 = jdut1 + jdut1f;
    
    % Convert from TEME to ECI
    [r_vectECI, v_vectECI, aeci] = teme2eci(r_vectTEME', v_vectTEME', 0, ttt, dpsi, deps);  % km, km/s
    
    % Keep the ECEF Position instead of ECI for the 3D Representation
    [r_vectECEF, v_vectECEF] = teme2ecef(r_vectTEME', v_vectTEME', 0, ttt, jdut1, lod, X_EOP, Y_EOP, 2);

    XECI_SGP4(k, :) = [r_vectECI' v_vectECI'];
    XECEF_SGP4(k, :) = [r_vectECEF' v_vectECEF'];

end

X0_SGP4 = XECI_SGP4(1, :)';
Xf_SGP4 = XECI_SGP4(end, :)';


%% Laplace Initial Orbit Determination

indices = [1 floor(N/2) N];      % choice of the observations set

% Initialize input State and Times
X0 = [];
times = [];
Parameters = [DAT, ut1_utc, X_EOP, Y_EOP, dX_EOP, dY_EOP];

% Define the Substate from the given indices
for i = 1 : length(indices)
    X0 = [X0; Z(i, :)];
    times = [times; tspan_obs(i)];
end

% We perform IOD at intermediate time -> t_star = times(2)
epoch = times(2);
[Rgs, dRgs, ddRgs] = lla2RdRddR(long_site, latgd_site, alt_site, epoch, Parameters);

% % Comparison with Older Method
% [Rgs_old, dRgs_old, ddRgs_old] = groundStation_Old(long_site, latgd_site, alt_site, epoch);
% fprintf(log, 'Rgs = [%.4f\t%.4f\t%.4f]\nRgs_old = [%.4f\t%.4f\t%.4f]\ndRgs = [%.4f\t%.4f\t%.4f]\ndRgs_old = [%.4f\t%.4f\t%.4f]\nddRgs = [%.4f\t%.4f\t%.4f]\nddRgs_old = [%.4f\t%.4f\t%.4f]\n', Rgs, Rgs_old, dRgs, dRgs_old, ddRgs, ddRgs_old);

[rLap, vLap] = LaplaceIOD(X0, times, Rgs, dRgs, ddRgs, epoch, mu);
coeLap = rvECI2coe(rLap, vLap, mu);

X0_LAP = [rLap; vLap];

% Display the Results
fprintf('\n\n\t<strong>Laplace Assumption Classical Orbital Elements:</strong>\n')
fprintf('\ta = %.4f\n\te = %.4f\n\ti = %.4f\n\tOmega = %.4f\n\tomega = %.4f\n\tnu0 = %.4f\n\n', [coeLap(1:2), rad2deg(coeLap(3:6))])
fprintf('\tThe reference altitude is: %.4f km\n', coeLap(1)-rE);

fprintf('\n\t<strong>State Comparison between Laplace Assumption and SGP4 Propagation:</strong>\n')
fprintf('\tX0_SGP4\t= [%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f]\n\tX0_LAP\t= [%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f]\n', XECI_SGP4(indices(2), :), X0_LAP);

if optional_plots
    figure('Name', 'Ground Station Location')
    ECIDrawTraj3D(Rgs'/norm(Rgs)*rE)
    if saveplots
        saveas(gcf, strcat('Output/Ground Station Location.jpg'))
    end
end


%% Laplace IOD with Target Inclination

if optimal_IOD

    % Define Sets of Observation indices
    all_indices = nchoosek(1:N, 3);     % choose any 3 out of N

    % Initialize Optimal COE
    coe_opt = [];
    target_diff = inf;

    % Initialize input State and Times
    X0 = [];
    times = [];
    P = [long_site, latgd_site, alt_site, mu]';

    % Try All Combinations
    for i = 1 : size(all_indices, 1)

        % Define the Substate from the given indices
        indices = all_indices(i, :);
        X0 = [];
        times = [];
        for j = 1 : length(indices)
            X0 = [X0; Z(indices(j), :)];
            times = [times; tspan_obs(indices(j))];
        end

        % We perform IOD at intermediate time -> t_star = times(2)
        epoch = times(2);
        [Rgs, dRgs, ddRgs] = lla2RdRddR(long_site, latgd_site, alt_site, epoch, Parameters);
        
        [rLap, vLap] = LaplaceIOD(X0, times, Rgs, dRgs, ddRgs, epoch, mu);
        coeLap = rvECI2coe(rLap, vLap, mu);

        % Update the Optimal Result
        if abs(coeLap(3) - target_incl) < target_diff
            target_diff = abs(coeLap(3) - target_incl);
            coe_opt = coeLap;
            indices_opt = indices;
        end

    end

    [X0_LAP(1:3), X0_LAP(4:6)] = coe2rvECI(coe_opt, mu);

    % Display the Optimization Results
    fprintf('\n\n\t<strong>Optimized Laplace Assumption with Target Set:</strong>\n')
    fprintf('\ta = %.4f\n\te = %.4f\n\ti = %.4f\n\tOmega = %.4f\n\tomega = %.4f\n\tnu0 = %.4f\n\n', [coe_opt(1:2), rad2deg(coe_opt(3:6))])
    fprintf('\tThe reference altitude is: %.4f km\n', coe_opt(1)-rE);
    fprintf('\tThe Optimal Order for the Observations is: %.0f %.0f %.0f \n', indices_opt)

    fprintf('\n\t<strong>State Comparison between Optimized Laplace Assumption and SGP4 Propagation:</strong>\n')
    fprintf('\tX0_SGP4\t= [%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f]\n\tX0_LAP\t= [%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f]\n', XECI_SGP4(indices_opt(2), :), X0_LAP);

end

%% Circular Hypothesis Initial Orbit Determination

indices = [1 floor(N/2)];        % choice of the observations set

% Initialize input State and Times
X0 = [];
times = [];

% Define the Substate from the given indices
for i = 1 : length(indices)
    X0 = [X0; Z(indices(i), :)];
    times = [times; tspan_obs(indices(i))];
end

Rgs1 = lla2eci([rad2deg(latgd_site), rad2deg(long_site), alt_site*1e3], datevec(times(1)), 'IAU-2000/2006', DAT, ut1_utc, [X_EOP, Y_EOP], 'dcip', [dX_EOP, dY_EOP])' * 1e-3;   % km
Rgs2 = lla2eci([rad2deg(latgd_site), rad2deg(long_site), alt_site*1e3], datevec(times(2)), 'IAU-2000/2006', DAT, ut1_utc, [X_EOP, Y_EOP], 'dcip', [dX_EOP, dY_EOP])' * 1e-3;   % km

[r1Circ, v1Circ] = CircularIOD(X0, times, Rgs1, Rgs2, mu, rho1_0);

coeCirc = rvECI2coe(r1Circ, v1Circ, mu);
X0_CIRC = [r1Circ; v1Circ];

% Display the Results
fprintf('\n\n\t<strong>Circular Assumption Classical Orbital Elements:</strong>\n')
fprintf('\ta = %.4f\n\te = %.4f\n\ti = %.4f\n\tOmega = %.4f\n\tomega = %.4f\n\tnu0 = %.4f\n\n', [coeCirc(1:2), rad2deg(coeCirc(3:6))])
fprintf('\tThe reference altitude is: %.4f km\n', coeCirc(1)-rE);

fprintf('\n\t<strong>State Comparison between Circular Assumption and SGP4 Propagation:</strong>\n')
fprintf('\tX0_SGP4\t= [%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f]\n\tX0_CIRC\t= [%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f]\n', XECI_SGP4(indices(1), :), X0_CIRC);


%% Circular IOD with Target Inclination

if optimal_IOD

    % Define Sets of Observation indices
    all_indices = nchoosek(1:N, 2);     % choose any 3 out of N

    % Initialize Optimal COE
    coe_opt = [];
    target_diff = inf;

    % Initialize input State and Times
    X0 = [];
    times = [];
    P = [long_site, latgd_site, alt_site, mu]';

    % Try All Combinations
    for i = 1 : size(all_indices, 1)

        % Define the Substate from the given indices
        indices = all_indices(i, :);
        X0 = [];
        times = [];
        for j = 1 : length(indices)
            X0 = [X0; Z(indices(j), :)];
            times = [times; tspan_obs(indices(j))];
        end

        Rgs1 = lla2eci([rad2deg(latgd_site), rad2deg(long_site), alt_site*1e3], datevec(times(1)), 'IAU-2000/2006', DAT, ut1_utc, [X_EOP, Y_EOP], 'dcip', [dX_EOP, dY_EOP])' * 1e-3;   % km
        Rgs2 = lla2eci([rad2deg(latgd_site), rad2deg(long_site), alt_site*1e3], datevec(times(2)), 'IAU-2000/2006', DAT, ut1_utc, [X_EOP, Y_EOP], 'dcip', [dX_EOP, dY_EOP])' * 1e-3;   % km
        
        [r1Circ, v1Circ] = CircularIOD(X0, times, Rgs1, Rgs2, mu, rho1_0);
        
        coeCirc = rvECI2coe(r1Circ, v1Circ, mu);

        % Update the Optimal Result
        if abs(coeCirc(3) - target_incl) < target_diff
            target_diff = abs(coeCirc(3) - target_incl);
            coe_opt = coeCirc;
            indices_opt = indices;
        end

    end

    [X0_CIRC(1:3), X0_CIRC(4:6)] = coe2rvECI(coe_opt, mu);

    % Display the Optimization Results
    fprintf('\n\n\t<strong>Optimized Circular Assumption with Target Set:</strong>\n')
    fprintf('\ta = %.4f\n\te = %.4f\n\ti = %.4f\n\tOmega = %.4f\n\tomega = %.4f\n\tnu0 = %.4f\n\n', [coe_opt(1:2), rad2deg(coe_opt(3:6))])
    fprintf('\tThe reference altitude is: %.4f km\n', coe_opt(1)-rE);
    fprintf('\tThe Optimal Order for the Observations is: %.0f %.0f %.0f \n', indices_opt)

    fprintf('\n\n\t<strong>State Comparison between Optimized Circular Assumption and SGP4 Propagation:</strong>\n')
    fprintf('\tX0_SGP4\t= [%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f]\n\tXf_KF\t= [%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f]\n', XECI_SGP4(indices_opt(1), :), X0_CIRC);

end


%% Kalman Filter from Circular Assumption

fprintf('\n\n\tApplying the Kalman Filter from Circular Assumption Results...\n')

% Define the Initial Covariance Matrix with TLE Accuracy
P0 = [sigma_r^2*eye(3), zeros(3, 3);...
      zeros(3, 3), sigma_v^2*eye(3)];

% Define Initial Conditions
X0 = X0_SGP4;

% Define Measurement Error Matrix
R = sigma_meas^2*eye(2);            % a.k.a. W_obs

% Define Model Error Matrix 
Q = 1e-6*eye(6);                    % a.k.a. W_apr
Q(1:3, 1:3) = zeros(3);

% Define a State Storage Matrix
X = zeros(6, N);

% Define a Covariance Storage Matrix
Ps = [];
Ps = [Ps, P0];
X(:,1) = X0;


% Kalmann Iterations
t0 = 0;
P_traces = zeros(N, 1);
P_traces = trace(P0);
P_est_traces = zeros(N, 1);

fprintf(log, '\n\nKalman Filter Iterations log:\n');

for j = 1 : N-1 

    if j > 1
        t0 = tf;
        X0 = Xf;
        P0 = Pk;
    end

    dt = seconds(tspan_obs(j+1) - tspan_obs(j));
    tf = t0 + dt;

    Phi0 = eye(6);
    X0 = [X0; reshape(Phi0, 36, 1)];

    % State Propagation
    [tspan_i, X_i] = ode113('DynamicalModel_Kalman', [t0, tf], X0, options);
    Xf = X_i(end, 1:6)';
    phi = X_i(end, 7:end)';
    Phi = reshape(phi, 6, 6);

    % Covariance Propagation
    P_est = Phi * P0 * Phi' + Q;
    P_est_traces(j+1) = trace(P_est);

    % Compute Kalman Gain Matrix
    Rgs = lla2eci([rad2deg(latgd_site), rad2deg(long_site), alt_site*1e3], datevec(tspan_obs(j+1)), 'IAU-2000/2006', DAT, ut1_utc, [X_EOP, Y_EOP], 'dcip', [dX_EOP, dY_EOP])' * 1e-3;       % km
    z = Z(j+1, :)';
    z_prev = Z(j, :)';
    H = H_matrix(Xf, Rgs);
    K = P_est * H' / ((H*P_est*H' + R));

    % Compute the Innovation and State Update

    % % Theoretically Correct Method - Not Working
    % Innovation = z - H*Xf;
    
    % Alternative Method 1
    [Raf, Decf] = ECI2TopoRaDec(Rgs, Xf(1:3));
    zf = [deg2rad(Raf) deg2rad(Decf)]';
    Innovation = z - zf;

    % % Alternative Method 2
    % Rgs0 = lla2eci([rad2deg(latgd_site), rad2deg(long_site), alt_site*1e3], datevec(tspan_obs(j)), 'IAU-2000/2006', DAT, ut1_utc, [X_EOP, Y_EOP], 'dcip', [dX_EOP, dY_EOP])' * 1e-3;       % km
    % [Ra0, Dec0] = ECI2TopoRaDec(Rgs0, X0(1:3));
    % z0 = [deg2rad(Ra0) deg2rad(Dec0)]';
    % Innovation = z - (z0 + H*(Xf-X0(1:6)));


    Xf_upd = K*Innovation;
    
    fprintf(log, 'Ra Innovation = %.6f\t\tDec Innovation = %.6f\n', Innovation);
    fprintf(log, 'Update to State Ratio: [%.3f %.3f %.3f %.3f %.3f %.3f]\n', Xf_upd./Xf);
    
    % Update of State and Covariance
    Xf = Xf + Xf_upd;
    Pk = (eye(6) - K*H) * P_est;

    % Debug - Check if Phi and P_est are positive definite
    if ~isPositiveDefinite(Pk) || ~isPositiveDefinite(P_est)
        error('Either Pk/k or Pk-1/k are NOT positive definite at iteration %d.', j);
    end

    % Store the Results
    P_traces(j+1) = trace(Pk);
    Ps = [Ps, Pk];

end


Xf_KF = Xf;


% Create the Visuals
P_vect = [];
its = [];
tspan_visuals = [];
for i = 1 : N

    if i == 1
        P_vect = [P_vect, P_traces(i)];
        its = [its, i];
        tspan_visuals = [tspan_visuals, tspan_obs(i)];
    end

    if i > 1 && i < N
        P_vect = [P_vect, P_est_traces(i), P_traces(i)];
        its = [its, i, i];
        tspan_visuals = [tspan_visuals, tspan_obs(i), tspan_obs(i)];
    end

    if i == N
        P_vect = [P_vect, P_traces(i)];
        its = [its, i];
        tspan_visuals = [tspan_visuals, tspan_obs(i)];
    end

end


COEKF = rvECI2coe(Xf(1:3), Xf(4:6), mu);
fprintf('\n\n\t<strong>Kalmann Filter Classical Orbital Elements:</strong>\n')
fprintf('\ta = %.6f\n\te = %.6f\n\ti = %.6f\n\tOmega = %.6f\n\tomega = %.6f\n\tnu0 = %.6f\n\n', [COEKF(1:2), rad2deg(COEKF(3:6))])
fprintf('\tThe reference altitude is: %.4f km\n', COEKF(1)-rE);

fprintf('\n\t<strong>State Comparison between Kalman Filter and SGP4 Propagation:</strong>\n')
fprintf('\tXf_SGP4\t= [%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f]\n\tXf_KF\t= [%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f]\n\n', Xf_SGP4, Xf_KF);


figure('name', 'Covariance Matrix Evolution')
plot(tspan_visuals, P_vect, 'LineStyle','-', 'Marker','.', 'LineWidth', 1.2, 'MarkerSize', 15, 'color', '#0d95fc', 'MarkerEdgeColor', '#590dfc');
title('tr$(P_k)$', 'Interpreter','latex', 'FontSize', 15)
if saveplots
    saveas(gcf, strcat('Output/Covariance Matrix Trace.jpg'))
end


if optional_plots
    figure('Name', 'Evolution of the Determinants of Estimated and Propagated Covariance Matrices')
    subplot(1, 2, 1)
    plot(tspan_obs, P_traces, 'LineStyle','-', 'Marker','.', 'LineWidth', 1.2, 'MarkerSize', 15, 'color', '#0d95fc', 'MarkerEdgeColor', '#590dfc');
    % plot(1:N, P_traces, 'LineStyle','-', 'Marker','.', 'LineWidth', 1.2, 'MarkerSize', 15, 'color', '#0d95fc', 'MarkerEdgeColor', '#590dfc');
    title('tr$(P_{k/k})$', 'Interpreter','latex', 'FontSize', 15)
    subplot(1, 2, 2)
    plot(tspan_obs, P_est_traces, 'LineStyle','-', 'Marker','.', 'LineWidth', 1.2, 'MarkerSize', 15, 'color', '#0d95fc', 'MarkerEdgeColor', '#590dfc');
    % plot(1:N, P_est_traces, 'LineStyle','-', 'Marker','.', 'LineWidth', 1.2, 'MarkerSize', 15, 'color', '#0d95fc', 'MarkerEdgeColor', '#590dfc');
    title('tr$(P_{k/k-1})$', 'Interpreter','latex', 'FontSize', 15)
    if saveplots
        saveas(gcf, strcat('Output/Estimated and Propagated Covariance Matrices Traces.jpg'))
    end
end


fclose(log);


%% Plot of the Results

% Laplace Assumption Results Propagation
COE0 = rvECI2coe(X0_LAP(1:3), X0_LAP(4:6), mu);

% Define the Time Domain
t0 = 0;
a0 = COE0(1);                % initial semi-major axis
T = 2*pi*sqrt(a0^3/mu);      % orbital period
tf = T;

% Perform the Propagation
[tspan, COE] = ode113('DynamicalModel_COE', [t0, tf], COE0, options);

% Store the Results
COE = COE(:, 1:6);
a = COE(:, 1);
e = COE(:, 2);
incl = COE(:, 3);
omega= COE(:, 4);
Omega = COE(:, 5);
nu = COE(:, 6);

M = length(tspan);
rMatrixECI_LAP = zeros(M, 3);

for j = 1 : M
    [r_vect, v_vect] = coe2rvECI(COE(j, :), mu);
    rMatrixECI_LAP(j, :) = r_vect';
end



% Circular Assumption Results Propagation
COE0 = rvECI2coe(X0_CIRC(1:3), X0_CIRC(4:6), mu);

% Define the Time Domain
t0 = 0;
a0 = COE0(1);                % initial semi-major axis
T = 2*pi*sqrt(a0^3/mu);      % orbital period
tf = T;

% Perform the Propagation
[tspan, COE] = ode113('DynamicalModel_COE', [t0, tf], COE0, options);

% Store the Results
COE = COE(:, 1:6);
a = COE(:, 1);
e = COE(:, 2);
incl = COE(:, 3);
omega= COE(:, 4);
Omega = COE(:, 5);
nu = COE(:, 6);

M = length(tspan);
rMatrixECI_CIRC = zeros(M, 3);

for j = 1 : M
    [r_vect, v_vect] = coe2rvECI(COE(j, :), mu);
    rMatrixECI_CIRC(j, :) = r_vect';
end



% Kalman Filter Results Propagation
COE0 = COEKF;

% Define the Time Domain
t0 = 0;
a0 = COE0(1);                % initial semi-major axis
T = 2*pi*sqrt(a0^3/mu);      % orbital period
tf = T;

% Perform the Propagation
[tspan, COE] = ode113('DynamicalModel_COE', [t0, tf], COE0, options);

% Store the Results
COE = COE(:, 1:6);
a = COE(:, 1);
e = COE(:, 2);
incl = COE(:, 3);
omega= COE(:, 4);
Omega = COE(:, 5);
nu = COE(:, 6);

tspan = tspan_obs(end) + seconds(tspan);

M = length(tspan);
rMatrixECI_KF = zeros(M, 3);

for j = 1 : M
    [r_vect, v_vect] = coe2rvECI(COE(j, :), mu);
    rMatrixECI_KF(j, :) = r_vect';
end


% SGP4 Results Propagation
COE0 = rvECI2coe(X0_SGP4(1:3), X0_SGP4(4:6), mu);
rMatrixECI_SGP4 = zeros(N, 3);

% Define the Time Domain
t0 = tspan_obs(1);
a0 = COE0(1);                % initial semi-major axis
T = 2*pi*sqrt(a0^3/mu);      % orbital period
tf = t0 + seconds(T);

% Compute tsince span
tsince0 = minutes(t0 - TLE_epoch);
tsincef = minutes(tf - TLE_epoch);
m = 1000;
tsince_span = linspace(tsince0, tsincef, m);


for k = 1 : m   % propagation in time

    % Compute tsince as time elapsed from TLE to Observation in minutes
    tsince = tsince_span(k);

    obs_epoch = t0 + minutes(tsince);
    obs_epoch_vect = datevec(obs_epoch);
    [obs_y, obs_mo, obs_d, obs_h, obs_mi, obs_s] = split_epochvec(obs_epoch_vect);

    % Propagation of position to Observation Time
    [~, r_vectTEME, v_vectTEME] = sgp4(satrec, tsince);
    
    % Evaluate ttt from UTC
    [~, ~, jdut1, jdut1f, ~, ~, ~, ttt, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~]...
        = convtime(obs_y, obs_mo, obs_d, obs_h, obs_mi, obs_s, 0, ut1_utc, DAT);     % timezone has been set to 0
    
    jdut1 = jdut1 + jdut1f;
    
    % Convert from TEME to ECI
    [r_vectECI, v_vectECI, aeci] = teme2eci(r_vectTEME', v_vectTEME', 0, ttt, dpsi, deps);  % km, km/s

    rMatrixECI_SGP4(k, :) = r_vectECI';

end


figure('Name', 'Trajectory Prediction in 3D for One Orbital Period ')
Lap_Traj = DrawTraj3D(rMatrixECI_LAP, '#22bf83');
Circ_Traj = DrawTraj3D(rMatrixECI_CIRC, '#227bbf');
KF_Traj = DrawTraj3D(rMatrixECI_KF);
SGP4_Traj = DrawTraj3D(rMatrixECI_SGP4, '#ff0505');

legend([Lap_Traj, Circ_Traj, KF_Traj, SGP4_Traj], {'Laplace Assumption', 'Circular Assumption', 'Kalman Filter Correction', 'SGP4 Propagation'}, 'Location', 'best', 'FontSize', 10)
if saveplots
    saveas(gcf, strcat('Output/Trajectories.jpg'))
end
