function Sats = ObsDuration(Sats)
% Description: this function fills the Obs_Duration variable in seconds of each
% Satellite with the time interval between first and last observations.


N = length(Sats);
M = length(Sats(1).Observable(:));

for i = 1 : N

    obs_initial = find(Sats(i).Observable(:), 1, "first");  % index of first observation

    streak = 0;
    
    for j = obs_initial+1 : M

        if Sats(i).Observable(j)

            streak = streak + 1;

        else

            break

        end

    end

    Sats(i).Obs_Duration = seconds(Sats(i).Epoch(obs_initial + streak) - Sats(i).Epoch(obs_initial));
    
    % obs_final = find(Sats(i).Observable(:), 1, "last");     % index of last observation

end




end