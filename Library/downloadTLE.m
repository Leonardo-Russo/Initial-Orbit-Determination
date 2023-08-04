function downloadTLE()
% Description: this function downloads the TLE Datasheet from the url
% defined locally.

load('Input/tle_data.mat');

if tle_date ~= datetime('today')
    tle_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle';
    tle_filename = 'Input/TLE.txt';
    websave(tle_filename, tle_url);
    tle_date = datetime('today');
    save('Input/tle_data.mat')
end

end