function Output = Concatenated_firing(Neurons, Ranges)

%This takes a Time ranges structure  from spike2, a Neurons structure containing
% all your spike times for each neuron,  and concatenates your firing across 
%all neurons in the time dimension. The output is a single vector
%containing all the times a spike occured during a state
% within the population of interest. 
%Written and contributed by Midha Ahmad, 12/12/2023.


% Inputs:
%            Neurons            A structure containing all your neurons and
%                               their properties
%            Ranges             Time ranges for state from Spike2

% Dependencies:
% Get_Events_In_Ranges.m
% Times_From_Ranges.m


x = [];
X = repmat(x,1,length(Neurons));
p = Times_From_Ranges(Ranges);
    for j = 1:length(Neurons);%loop through each neuron
        t = Neurons(j).times;
        events_in_ranges = Get_Events_In_Ranges(t,p);
        Firing_State = t(logical(events_in_ranges));
        X{j} = [Firing_State]; %firing rates for this state
        Xsz(j,:) = size(X{j}); % returns size of each vector
    end;

Colmax = max(Xsz(:,1)); % Maximum # of Columns
Statemtx = NaN(j,Colmax);% Preallocate/initalize
    for k1 = 1:j;
       Statemtx(k1,1:Xsz(k1,1)) = X{k1}.';% Fill Matrix
    end

Statemtx(isnan(Statemtx))=0; %replace nan values with zero in the summed matrix

%sort all the firings
State_FR_linear =  Statemtx(:);
State_FR_linear = sort(State_FR_linear);
%eliminate zero cells
State_FR_linear = nonzeros(State_FR_linear);
Output = State_FR_linear; %returns the population spiking  
end