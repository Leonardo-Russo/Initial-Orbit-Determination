function rad_angle = arcsec2rad(arcsec_angle)
% Description: this function converts the input angle from arcsec to
% radians.

rad_angle = arcsec_angle * pi / (180*3600);

end