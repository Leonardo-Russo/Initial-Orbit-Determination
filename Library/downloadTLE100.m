function downloadTLE100()
% Description: this function downloads the TLE Datasheet from the url
% defined locally.

load('Input/tle_data.mat');

if tle100_date ~= datetime('today')
    tle100_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=visual&FORMAT=tle';
    tle100_filename = 'Input/TLE100.txt';
    websave(tle100_filename, tle100_url);
    tle100_date = datetime('today');
    save('Input/tle_data.mat')
end

end