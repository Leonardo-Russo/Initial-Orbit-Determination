function RaH = Ra2HHMMSS(Ra)

str = char(duration((Ra / 360) * 24, 0, 0, 0));

RaH = [str2num(str(1:2)), str2num(str(4:5)), str2num(str(7:8))];

end