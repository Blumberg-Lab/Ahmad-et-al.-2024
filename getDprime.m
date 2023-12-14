function [dprime accuracy nSamples] = getDprime(setNumbers);

tP = setNumbers(1); fN = setNumbers(2);
fP = setNumbers(3); tN = setNumbers(4);
nSamples = sum(setNumbers);

hRate = tP / (tP + tN);
mRate = 1 - hRate;

faRate = fP / (fP + tN);
crRate = 1 - faRate;

accuracy = (tP + tN) / nSamples;
z_hit = norminv(hRate);
z_fa = norminv(faRate);

dprime = z_hit - z_fa;