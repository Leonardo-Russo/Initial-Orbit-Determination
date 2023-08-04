function [EOP] = get_EOP()

eopfile = 'EOP-Last5Years.txt';
EOPData = readlines(eopfile);

EOPstr = [];
loading_flag = 0;

for j = 1 : length(EOPData)
    
    line = EOPData(j);
    
    % Note that we do NOT distinguish between Observed and Predicted EOP's
    if line == "END OBSERVED" || line == "END PREDICTED"
        loading_flag = 0;
    end

    if loading_flag == 1
        EOPstr = [EOPstr; line];
    end

    if line == "BEGIN OBSERVED" || line == "BEGIN PREDICTED"
        loading_flag = 1;
    end

end

% Build EOP Matrix
EOP = [];

for j = 1 : length(EOPstr)
    eop_line = split(EOPstr(j))';
    eop_line = str2double(eop_line);
    EOP = [EOP; eop_line];
end

end