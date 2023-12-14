function [Output] = Firing_Rate_Violinplots(States, Neurons, norm)
% The following function is a multi-purpose visualization tool. You can visualize:
% 1. Average firing rate per neuron in each state as bar graphs.
% 2. Violin plots of avg firing rate, all neurons pooled together, with state as categories
% 2. Normalized versions (to max firing across states) of both plots.
% 
% Inputs:
%           States    a structure containing TIMES of of all states, imported from spike2  
%           Neurons   a structure containing all of neurons 
%           norm      1 = normalize each neuron to its max firing rate
%                         across states
%                     0 = do not normalize
% Outputs:
%              a structure containing mean FRs and normalized mean FRs for
%              each neuron per state
% Dependencies:
%              this code relies on the violinplot.m function by Bastian Bechtold, 2016. the 
%              script can be downloaded here:
%              https://github.com/bastibe/Violinplot-Matlab/blob/master/violinplot.m
% Parameters:
%              see parameters for the violins, whisker plots, scatter etc. 
%              in the violinplot.m function. They can be manually adjusted
%              below in the last section.
%
% Written and contributed by Midha Ahmad, 12/14/2023.





%grab the fieldnames of the states located in the structure
fnStates = fieldnames(States);

%first grab each neuron's firing per state
for k =1 : numel(Neurons)
    for s =1 : numel(fnStates) 
    eventRate.(string(Neurons(k).title)).(string(fnStates(s))) = Event_Rate(...
    States.(string(fnStates(s))), Neurons(k).times);
    end
end

%grabs each neuron's maximal firing rate across states
fnBar = fieldnames(eventRate);
for k = 1 : numel(fnBar)
    for s =1 : numel(fnStates)
     stderror.(string(fnBar(k)))(s) = std(eventRate.(string(fnBar(k))).(string(fnStates(s))).rates)...
     / sqrt(length(eventRate.(string(fnBar(k))).(string(fnStates(s))).rates));
     %grabs stderrors as well for plotting
     eventRate.(string(fnBar(k))).(string(fnStates(s))).rates = eventRate.(string(fnBar(k))).(string(fnStates(s))).rates;
     meanFRs.(string(fnBar(k)))(s) = nanmean(eventRate.(string(fnBar(k))).(string(fnStates(s))).rates);
     %this just takes the maximal firing for each state
     eventRate.(string(fnBar(k))).(string(fnStates(s))).maxFR = max(eventRate.(string(fnBar(k))).(string(fnStates(s))).rates);
    end
end



if norm == 1;
    for k = 1 : numel(fnBar)
        for s =1 : numel(fnStates)
         maxfirings.(string(fnBar(k)))(s) = eventRate.(string(fnBar(k))).(string(fnStates(s))).maxFR;
         m.(string(fnBar(k))) = max(maxfirings.(string(fnBar(k)))(s))
         %this takes the highest firing of neuron across entire run
         %now normalize that for the value you have grabbed above
         eventRate.(string(fnBar(k))).(string(fnStates(s))).normFR = eventRate.(string(fnBar(k))).(string(fnStates(s))).rates/m.(string(fnBar(k)));
         meannormFR.(string(fnBar(k)))(s) = nanmean(eventRate.(string(fnBar(k))).(string(fnStates(s))).normFR)
        end
    end
end

%% As a sanity check, this function generates each neuron's firing rate
%%bar graphs to make sure there is nothing alarming in indiviudal data


stateNames = categorical(fnStates);

if norm == 0;
    %plot bar graphs for all the neurons, no normalization
    figure
    title('Averaged firing rate across states')
    for h = 1 : length(fnBar)
    
        subplot(5,7,h)
        bar(stateNames,meanFRs.(string(fnBar(h))), 'BarWidth', 0.3)    
        
        hold on 
        
        er = errorbar(stateNames,meanFRs.(string(fnBar(h))),stderror.(string(fnBar(h))));    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';
        
        xlabel('States')
        ylabel('Firing Rate')
        title([fnBar(h)])
    end
end;

if norm == 1;
    %plot the normalized figure
    figure
    title('normalized avg firing rate across states')
    for h = 1 : length(fnBar)
        subplot(5,7,h)
        bar(stateNames,meannormFR.(string(fnBar(h))), 'BarWidth', 0.3)
        xlabel('States')
        ylabel('Firing Rate')
        title([fnBar(h)])
        ylim([0 1])
    end
end

%% Generate the violin plots%%


%convert your mean firing rate structure into a readable data
Data_violin = zeros(length(Neurons),3);
norm_Data_violin = zeros(length(Neurons),3);
cat_violin = repmat((stateNames'), [length(Neurons),1]);
    for h = 1 : length(fnBar)
        Data_violin(h,:) = meanFRs.(string(fnBar(h)));
        if norm == 1;
        norm_Data_violin(h,:) = meannormFR.(string(fnBar(h)));
        end
    end


%show show show
if norm == 0;
    figure
    violinplot(Data_violin, cat_violin,...
        'MarkerSize', 10,...,
        'ShowMean', false, ...
        'ShowMedian', true,...
        'ShowWhiskers', true,...
        'ShowBox', false,....
        'ShowNotches',false,...
        'ShowData', true,...
        'DataStyle', 'scatter')
    xlabel('States')
    ylabel('Firing Rate')
    title('Averaged firing rate of PZ neurons across states')

else
    figure
    violinplot(norm_Data_violin, cat_violin,...
        'MarkerSize', 10,...,
        'ShowMean', false, ...
        'ShowMedian', true,...
        'ShowWhiskers', true,...
        'ShowBox', false,....
        'ShowNotches',false,...
        'ShowData', true,...
        'DataStyle', 'scatter')
    xlabel('States')
    ylabel('Firing Rate')
    title('normalized firing rate of PZ neurons across states')
end

if norm == 0;
Output.meanFrs = meanFRs;
else
Output.normFrs = meannormFR;
end
end