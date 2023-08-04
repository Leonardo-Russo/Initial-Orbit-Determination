%% Beidou - Data Builder

close all
clear all
clc

filename = 'TDM_BEIDOU_38091_2022_11_02_SCUDO.kvn';
lines = readlines(filename);
l = length(lines);

tbias = -0.5756;    % s

read = 0;
data = [];

for i = 1 : l

    if lines(i) == "DATA_STOP"
        read = 0;
    end

    if read == 1
        data = [data; lines(i)];
    end

    if lines(i) == "DATA_START"
        read = 1;
    end

end

n = length(data);
RA = [];
DEC = [];
tspan_obs = [];

for i = 1 : n
    
    splits = strsplit(data(i));

    if splits(1) == "ANGLE_1"
        RA = [RA; deg2rad(double(splits(4)))];
        date = strsplit(splits(3), 'T');
        date = strcat(date(1), {' '}, date(2));
        time = datetime(date, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
        tspan_obs = [tspan_obs; time + seconds(tbias)];
    end

    if splits(1) == "ANGLE_2"
        DEC = [DEC; deg2rad(double(splits(4)))];
    end

end

m = n/2;
Z = zeros(m, 2);

for j = 1 : m
    Z(j, :) = [RA(j), DEC(j)];
end

N = length(tspan_obs);


save('Beidou.mat')