function EOPf = findEOPf(tf)
% Description: This function evaluates the EOP at a specified observation
% time provided as input.

% Reading EOP
EOP = get_EOP();

year = tf(1);
month = tf(2);
day = tf(3);

n = length(EOP);

for i = 1 : n

    if EOP(i, 1) == year && EOP(i, 2) == month && EOP(i, 3) == day
        EOPf = EOP(i, :);
        break
    end

end

EOPf = EOPf';

end