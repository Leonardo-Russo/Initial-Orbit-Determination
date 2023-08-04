%% Sentinel 3B - Data Builder

close all
clear all
clc

% Optical Measurements Sentinel 3B
RA1 = deg2rad(220.2384883);     % rad
DEC1 = deg2rad(-22.2310199);    % rad
RA2 = deg2rad(216.0834569);     % rad
DEC2 = deg2rad(-14.9416624);    % rad
RA3 = deg2rad(210.8583380);     % rad
DEC3 = deg2rad(-5.6574439);     % rad
RA4 = deg2rad(204.2886114);     % rad
DEC4 = deg2rad(5.6756058);      % rad
RA5 = deg2rad(195.8338798);     % rad
DEC5 = deg2rad(18.6624492);     % rad
RA6 = deg2rad(185.0851246);     % rad
DEC6 = deg2rad(31.5363487);     % rad
RA7 = deg2rad(171.9800163);     % rad
DEC7 = deg2rad(42.1061675);     % rad
RA8 = deg2rad(157.2297416);     % rad
DEC8 = deg2rad(49.1704531);     % rad
RA9 = deg2rad(142.6389246);     % rad
DEC9 = deg2rad(52.8259866);     % rad


tbias = -0.2479;    % s
texp = 1;           % s

tbias = tbias + texp/2;

% Times of Measurements (UTC)
t1 = datetime(2022, 06, 21, 21, 04, 41, 180+tbias);
t2 = datetime(2022, 06, 21, 21, 05, 20, 649+tbias);
t3 = datetime(2022, 06, 21, 21, 06, 00, 716+tbias);
t4 = datetime(2022, 06, 21, 21, 06, 40, 591+tbias);
t5 = datetime(2022, 06, 21, 21, 07, 20, 775+tbias);
t6 = datetime(2022, 06, 21, 21, 08, 00, 957+tbias);
t7 = datetime(2022, 06, 21, 21, 08, 41, 024+tbias);
t8 = datetime(2022, 06, 21, 21, 09, 21, 083+tbias);
t9 = datetime(2022, 06, 21, 21, 10, 00, 651+tbias);

% Define the Observations Matrix
tspan_obs = [t1 t2 t3 t4 t5 t6 t7 t8 t9]';
Z = [[RA1 RA2 RA3 RA4 RA5 RA6 RA7 RA8 RA9]',...
     [DEC1 DEC2 DEC3 DEC4 DEC5 DEC6 DEC7 DEC8 DEC9]'];

% n° of Measurements
N = length(tspan_obs);

save('Sentinel3B.mat')