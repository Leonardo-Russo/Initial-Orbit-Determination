function [tle1, tle2] = get_TLE(tlefile)
% Description: this function is AD HOC for this application, because
% normally the TLE also contain the name string, while in this case this
% line is NOT present!

TLEData = readlines(tlefile);   % read each line of input file

tle1 = char(TLEData(1));
tle2 = char(TLEData(2));


end