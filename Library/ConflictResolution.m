function Sats = ConflictResolution(Sats, dt)
% Description: this function returns a Sats Data Structure in which all
% conflicts have been solved.
% Planning Problem
N = length(Sats);
M = length(Sats(1).Epoch(:));

SatsGEO = [];
SatsLEO = [];

countersGEO = [];
indicesGEO = [];
countersLEO = [];
indicesLEO = [];

if dt == 60
    GEO_threshold = 6;
    LEO_threshold = 6;
elseif dt == 1
    GEO_threshold = 5*60;
    LEO_threshold = 5*60;
else
    error("Please, use either a time step of 1s or 60s.\n")
end

GEO_enough = 2*GEO_threshold;
LEO_enough = 2*LEO_threshold;


% Compute Local Variables
for i = 1 : N
    
    Sat = Sats(i);

    if Sat.SemiMajorAxis > 2

        SatsGEO = [SatsGEO; Sat];   % store satellite in GEOs

        index = find(Sat.Observable(:), 1, "first");   % index of first observation
        obs_counter = sum(Sat.Observable(:));
        
        countersGEO = [countersGEO; obs_counter];
        indicesGEO = [indicesGEO; index];
        
    else

        SatsLEO = [SatsLEO; Sat];   % repeat process for LEOs

        index = find(Sat.Observable(:), 1, "first");
        obs_counter = sum(Sat.Observable(:));
    
        countersLEO = [countersLEO; obs_counter];
        indicesLEO = [indicesLEO; index];

    end

end


NGEO = length(SatsGEO);
NLEO = length(SatsLEO);


% Remove GEO Sats under Observability Threshold
i = 1;
while true

    clc
    fprintf("Applying Observability Threshold to GEO: %.0f / %.0f\n", [i, NGEO])
    
    if countersGEO(i) < GEO_threshold
        SatsGEO = RemoveItem(SatsGEO, i);
        countersGEO = RemoveItem(countersGEO, i);
        indicesGEO = RemoveItem(indicesGEO, i);
        i = i-1;
        NGEO = length(SatsGEO);
    end

    if i >= NGEO
        break
    end
    
    i = i+1;

end

% Remove LEO Sats under Observability Threshold
i = 1;
while true
    
    clc
    fprintf("Applying Observability Threshold to LEO: %.0f / %.0f\n", [i, NLEO])

    if countersLEO(i) < LEO_threshold
        SatsLEO = RemoveItem(SatsLEO, i);
        countersLEO = RemoveItem(countersLEO, i);
        indicesLEO = RemoveItem(indicesLEO, i);
        i = i-1;
        NLEO = length(SatsLEO);
    end

    if i >= NLEO
        break
    end
    
    i = i+1;

end


NGEO = length(SatsGEO);
NLEO = length(SatsLEO);

observedGEO = zeros(NGEO, 1);
observedLEO = zeros(NLEO, 1);


%%% Conflict Resolution for GEO %%%
for j = 1 : M
    
    clc
    fprintf("Conflict Resolution for GEOs at epoch %.0f / %.0f\n", [j, M])

    % Compute the Contenders
    contenders = [];
    lengths = [];
    indices = [];
    for i = 1 : NGEO
        if SatsGEO(i).Observable(j) == 1
            contenders = [contenders; SatsGEO(i)];
            lengths = [lengths; countersGEO(i)];     % subset of counter for contenders
            indices = [indices; i];
        end
    end
    
    C = length(contenders);

    [~, lengths_sortidx] = sort(lengths, 'ascend');

    for k = 1 : C   % finds the most worthy
        candidate = lengths_sortidx(k);
        if observedGEO(candidate) < GEO_enough
            win_index = candidate;
            break
        end

        if k == C
            win_index = lengths_sortidx(1);
            break
        end
    end

    winner = contenders(win_index);
    
    losers_count = 0;
    for r = 1 : C
        if r ~= win_index
            SatsGEO(indices(r)).Observable(j) = 0;
            losers_count = losers_count + 1;
        else
            observedGEO(indices(r)) = observedGEO(indices(r)) + 1;
        end
    end
    
end

SatsGEO = RemoveUnobservables(SatsGEO);


%%% Conflict Resolution for LEO %%%
for j = 1 : M
    
    clc
    fprintf("Conflict Resolution for LEOs at epoch %.0f / %.0f\n", [j, M])

    % Compute the Contenders
    contenders = [];
    lengths = [];
    indices = [];
    streaks = [];
    for i = 1 : NLEO
        if SatsLEO(i).Observable(j) == 1
            contenders = [contenders; SatsLEO(i)];
            lengths = [lengths; countersLEO(i)];     % subset of counter for contenders
            indices = [indices; i];
            streaks = [streaks; LEO_threshold - observedLEO(i)];
        end
    end
    
    C = length(contenders);
    
    [~, lengths_sortidx] = sort(abs(lengths-mean(lengths)) + 5*streaks, 'ascend');
    

    for k = 1 : C   % finds the most worthy
        candidate = lengths_sortidx(k);
        if observedLEO(candidate) < LEO_enough
            win_index = candidate;
            break
        end

        if k == C
            win_index = lengths_sortidx(1);
            break
        end
    end

    winner = contenders(win_index);
    
    losers_count = 0;
    for r = 1 : C
        if r ~= win_index
            SatsLEO(indices(r)).Observable(j) = 0;
            losers_count = losers_count + 1;
        else
            observedLEO(indices(r)) = observedLEO(indices(r)) + 1;
        end
    end
    
end

SatsLEO = RemoveUnobservables(SatsLEO);


NGEO = length(SatsGEO);
NLEO = length(SatsLEO);


% Join LEOs and GEOs Satellites Prioritizing LEOs
for j = 1 : M

    for iL = 1 : NLEO
        free = 1;
        if SatsLEO(iL).Observable(j) == 1
            free = 0;
            break
        end
    end

    if ~ free
        for iG = 1 : NGEO
            if SatsGEO(iG).Observable(j) == 1
                SatsGEO(iG).Observable(j) = 0;
            end
        end
    end

end


SatsGEO = RemoveUnobservables(SatsGEO);

Sats = [SatsGEO; SatsLEO];


end