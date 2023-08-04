%% Sentinel 3A - Data Builder

close all
clear all
clc

% Optical Measurements Sentinel 3A
RA1 = deg2rad(210.5091647);     % rad
DEC1 = deg2rad(-18.0950098);    % rad
RA2 = deg2rad(205.6703196);     % rad
DEC2 = deg2rad(-10.4010653);    % rad
RA3 = deg2rad(199.8210477);     % rad
DEC3 = deg2rad(-1.2152585);     % rad
RA4 = deg2rad(192.8592320);     % rad
DEC4 = deg2rad(9.0517257);      % rad
RA5 = deg2rad(184.4013246);     % rad
DEC5 = deg2rad(19.8849237);     % rad
RA6 = deg2rad(174.6870341);     % rad
DEC6 = deg2rad(29.7102390);     % rad
RA7 = deg2rad(163.9949752);     % rad
DEC7 = deg2rad(37.5376987);     % rad

tbias = -0.3136;    % s
texp = 1;           % s

tbias = tbias + texp/2;

% Times of Measurements (UTC)
t1 = datetime(2022, 06, 22, 21, 18, 01, 062+tbias);
t2 = datetime(2022, 06, 22, 21, 18, 40, 723+tbias);
t3 = datetime(2022, 06, 22, 21, 19, 20, 728+tbias);
t4 = datetime(2022, 06, 22, 21, 20, 00, 412+tbias);
t5 = datetime(2022, 06, 22, 21, 20, 40, 856+tbias);
t6 = datetime(2022, 06, 22, 21, 21, 20, 789+tbias);
t7 = datetime(2022, 06, 22, 21, 22, 00, 350+tbias);

% Define the Observations Matrix
tspan_obs = [t1 t2 t3 t4 t5 t6 t7]';
Z = [[RA1 RA2 RA3 RA4 RA5 RA6 RA7]',...
     [DEC1 DEC2 DEC3 DEC4 DEC5 DEC6 DEC7]'];

% n° of Measurements
N = length(tspan_obs);

save('Sentinel3A.mat')