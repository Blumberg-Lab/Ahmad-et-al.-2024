function ER = Event_Rate(ranges,events)
%The following function takes time ranges, and calcualtes the rate at which
%an event occured in those times.

% Inputs:
%            events            All the times an event happened occured across
%                              your run
%            Ranges            Time ranges for state from Spike2

%Dependencies
% NONE

%Contributed by Midha Ahmad, 12/13/2023. Adapted from Kev You.

%% Parameters
Params = inputParser;
Params.addRequired('ranges',@(x) isnumeric(x) || islogical(x));
Params.addRequired('events',@(x) isnumeric(x) || islogical(x));


%make new array with all the periods of interest
period_idx = cell(1, 1000);

for i = 1 : 2 : length(ranges)
    
        period_idx{i} = [ranges(i), ranges(i+1)];
        
        if ranges(i+1) - ranges(i) < 0
            
            period_idx{i} = {};
        end

end

%eliminate empty cells
period_idx = period_idx(~cellfun(@isempty, period_idx));
event_times = []

%capture event times
for n = 1 : length(period_idx)
    eventCount = 0;
for m = 1 : length(events)
    currentCell = period_idx(n);
    if events(m) >= currentCell{1}(1) && events(m) <= currentCell{1}(2)
        eventCount = eventCount + 1;
            event_times{n} = [events(m)];
    else 
        eventCount = eventCount;
    end
end
   eventRate(n) = eventCount./(currentCell{1}(2) - currentCell{1}(1)); %calculates the rate

end


%outputs in a structure format
ER.rates = eventRate;
ER.periods = period_idx;
ER.times = event_times;
ER.eventcount = eventCount;
end
