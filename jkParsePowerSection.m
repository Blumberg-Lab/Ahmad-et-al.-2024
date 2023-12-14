function [datSigLogical datSigTS] = jkParsePowerSection(inputData, inputDataTS, SF, threshold, minDuration, smoothingFlag)
%-----ENDOF

if size(inputData, 1) < size(inputData, 2)
    inputData = transpose(inputData);
end %size(inputData, 1) < size(inputData, 2)

envData = envelope(inputData, 500, 'analytic');
threshold = median(envData) * threshold;
if smoothingFlag == 0
    byThreshold = envData > threshold;
else
    byThreshold = smooth(envData, 500) > threshold;
end %smoothingFlag == 0
byDuration = zeros(length(byThreshold), 1);
byThresholddiff = [0 ; diff(byThreshold)]; stIdx2BD = find(byThresholddiff == 1); edIdx2BD = find(byThresholddiff == -1);
for stRUN = 1:1:length(stIdx2BD)
    thisST = stIdx2BD(stRUN); thisED = edIdx2BD(min(find(edIdx2BD > thisST))) - 1;
    if (thisST < thisED) & all(byThreshold(thisST:thisED) == 1) & ((thisED - thisST + 1) > minDuration * SF)
        byDuration(thisST:thisED) = 1;
    end %(thisST < thisED) & all(byThreshold(thisST:thisED) == 1) & ((thisED - thisST + 1) > minDuration * SF)
end %stRUN = 1:1:length(stIdx2BD)
datSigLogical = (byThreshold & byDuration);

datSigLogicalDiff = [0 ; diff(datSigLogical)]; stIdx2BD = find(datSigLogicalDiff == 1); edIdx2BD = find(datSigLogicalDiff == -1);
datSigTS = [];
for stRUN = 1:1:length(stIdx2BD)
    thisST = stIdx2BD(stRUN); thisED = edIdx2BD(min(find(edIdx2BD > thisST))) - 1;
    
    if (thisST < thisED) & all(datSigLogical(thisST:thisED) == 1)
        datSigTS = [datSigTS ; ... 
                    [inputDataTS(thisST) inputDataTS(thisED)]];
    end %(thisST < thisED) & all(datSigLogical(thisST:thisED) == 1)
end %stRUN = 1:1:length(stIdx2BD)