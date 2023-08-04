function angle = checkAngle(angle, type)
% Description: this function brings any input angle into the 0°-360°
% interval.

if type == 'rad'
    ref = 2*pi;
elseif type == 'deg'
    ref = 360;
else
    error('Angle Unit not Valid!')
end

while true
        
    if angle < 0
        angle = angle + ref;
    elseif angle > ref
        angle = angle - ref;
    else
        break
    end

end


end