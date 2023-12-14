function [midAligned] = jkGetMidAliSignals(givenSignal, thisMidPoint, nBinsPrePost)
%
%[givenSignal]
%[thisMidPoint]
%[nBinsPrePost] how many samples will be put before and after the mid-point [if data were shorter, Nans will be filled in instead]
%
%[midAligned] will be a [1 x (nBinsPrePost*2 + 1)] matrix
%
%June-18-2023, Jangjin Kim

midAligned = [];

thisFullLength = length(givenSignal);
    if size(givenSignal, 1) > size(givenSignal, 2)
        givenSignal = transpose(givenSignal);
    end %size(givenSignal, 1) > size(givenSignal, 2)

%pre-
thisPre = [];
if thisMidPoint > nBinsPrePost
    thisPre = givenSignal((thisMidPoint-nBinsPrePost):thisMidPoint);
else
    nNans = 1 + (nBinsPrePost - thisMidPoint);
    thisPre = [nan(1, nNans) givenSignal(1:thisMidPoint)];
end %thisMidPoint > nBinsPrePost

%post-
thisPost = [];
if thisFullLength >= (thisMidPoint + nBinsPrePost - 1)
    thisPost = givenSignal(thisMidPoint+1:thisMidPoint+nBinsPrePost-1);
else
    nNans = (thisMidPoint + nBinsPrePost - 1) - thisFullLength;
    thisPost = [givenSignal(thisMidPoint+1:end) nan(1, nNans)];
end %thisFullLength >= (thisMidPoint + nBinsPrePost - 1)
    
midAligned = [thisPre thisPost];