function [Output] = Silence_bouts(State_FR_linear, mincrit,maxcrit,bin)

%This function uses all the times a spike occured within the population, to
%calculate the silences that existed at a population-level. NOTE: This function
%can also be adapted to caluclate inter-spike intervals, depending on the 
%critera you define.

%for sanity checks, this function creates two figures:
%     1. a figure of all the state-dependent spiking within a population
%        across the whole run.     
%     2. a frequency histogram distribution of the intervals output,
%        normalized by probablity 


% Inputs:
%            State_FR_linear    Output vector from Concatenated_Firing.m
%            mincrit            time interval, in seconds, below which any 
%                               durations will not be considered
%            maxcrit            time interval, in seconds, above which any 
%                               durations will not be considered
%            bin                bin size, in seconds, of the histogram
% Output:
%           a structure containing silence-bout histogram values, and the silence
%           bout durations.

% Dependencies:
% Get_Events_In_Ranges.m
% Times_From_Ranges.m
% Concatenated_Firing.m

%Written and contributed by Midha Ahmad, 12/13/2023.


figure
scatter(State_FR_linear, zeros(size(State_FR_linear)),100,'Color',[0 0 1], 'Marker', "|")
% note, above plots all spikes in one command

%calculate inter-bout/inter-spike intervals
interspikes_all = diff(State_FR_linear);



%bin size of all the silences greater than mincrit
longintervals = interspikes_all > mincrit; %gets rid of interrupted bouts
longintervals2 = interspikes_all(longintervals);
longintervals = longintervals2 < maxcrit;%gets rid of long silences
longintervals3 = longintervals2(longintervals);

Stateintervals = longintervals3;

%make a frequency histogram distribution of the intervals
dt = bin; % in s, because spike times are in s
isi_edges = 0:dt:1; % bin edges for ISI histogram
 
figure
[h] = histogram(Stateintervals,isi_edges,Normalization="probability")
xlim([mincrit 1])
xlabel('ISI of State silences')
HistogramVal = h.Values;
ax = h.BinEdges;


Output.Hval = HistogramVal; %y-axis
Output.freqbins = ax;%x-axis
Output.intervals = Stateintervals;%all the durations that met your criteria

end