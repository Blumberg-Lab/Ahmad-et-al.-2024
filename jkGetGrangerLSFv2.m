function [GCMat time2use freq2use timeOrder bicMat] = jkGetGrangerLSFv2(sig1, sig2, FreqRange, FreqScale)
%[GCMat time2alz freq2use timeOrder] = jkComputeGranger(sig1, sig2, SF, FreqRange)
%Algorithms were adopted from bsmart toollbox [make sure to refer the article]
%
%Code by Jangjin Kim, Sep-18-2020

%Check the validity of params
if ~(size(sig1, 1) == size(sig2, 1) & size(sig1, 2) == size(sig2, 2))
	error(['signal#1 and #2 should be the same size']);
end	%~(size(sig1, 1) == size(sig2, 1) & size(sig1, 2) == size(sig2, 2))

%SF
SF = 30.5180;

%Define time space
stEOI = 28;
edEOI = 137;

szEOI = 55;
thisBuff = 27;

timeSpace = linspace(-1000, 5000, size(sig1, 2));
time2use = [0 4000];
time2useIdx = [stEOI edEOI];

time2alz = time2useIdx(1):thisBuff:time2useIdx(2);

%Define freq space
nFreqSpace = 20;		%35 in previous versions
if FreqScale == 1
	freq2use = linspace(min(FreqRange), max(FreqRange), nFreqSpace);
else
	freq2use = logspace(log10(min(FreqRange)), log10(max(FreqRange)), nFreqSpace);
end %FreqScale == 1

nr = size(sig1, 1);
nl = szEOI;

GCMat = zeros(length(freq2use), length(time2alz), 2);		%3rd dim for x->y and y->x

%Remove ERP component
sig1 = bsxfun(@minus, sig1, nanmean(sig1, 1));
sig2 = bsxfun(@minus, sig2, nanmean(sig2, 1));

minBIC = 3; maxBIC = 12; nBicSpace = 6; bicMat = nan(length(time2alz), nBicSpace);
if maxBIC <= minBIC
	error(['EOI is too small']);
else
	bicSpace = linspace(minBIC, maxBIC, nBicSpace);
end %maxBIC <= minBIC

for timeRUN = 1:1:length(time2alz)
    dat1_t = sig1(:, (time2alz(timeRUN)-thisBuff):(time2alz(timeRUN)+thisBuff));
    dat2_t = sig2(:, (time2alz(timeRUN)-thisBuff):(time2alz(timeRUN)+thisBuff));

	%detrending and zscoring
	for trialRUN = 1:1:size(dat1_t, 1)
		dat1_t(trialRUN, :) = zscore(detrend(dat1_t(trialRUN, :)));
		dat2_t(trialRUN, :) = zscore(detrend(dat2_t(trialRUN, :)));
	end %trialRUN = 1:1:size(dat1_t, 1)

	dat1_t = reshape(transpose(dat1_t), 1, nr * nl);
	dat2_t = reshape(transpose(dat2_t), 1, nr * nl);

	dat_t = [dat1_t ; dat2_t];

	%Get BIC
	for bicRUN = 1:1:nBicSpace
		[Axy E] = armorf(dat_t, nr, nl, bicSpace(bicRUN));
		bicMat(timeRUN, bicRUN) = log(det(E)) + (log(nr * nl) * bicSpace(bicRUN) * 2 ^2) / (nr * nl);

		clear Axy E;
	end %bicRUN = 1:1:nBicSpace

	clear dat1_t dat2_t dat_t;
end %timeRUN = 1:1:length(time2alz)

%Detect the point where the estimation was significantly improved
avgBicMat = nanmean(bicMat, 1);
BICpercentChange = diff(avgBicMat) ./ nanmean(bicMat(:, 1), 1);
[dummy idxMaxChange] = nanmax(BICpercentChange);

szOrder = bicSpace(1, (1 + idxMaxChange));
timeOrder = szOrder * (1000 / SF);

%Get Granger with optimized BIC
for timeRUN = 1:1:length(time2alz)
    dat1_t = sig1(:, (time2alz(timeRUN)-thisBuff):(time2alz(timeRUN)+thisBuff));
    dat2_t = sig2(:, (time2alz(timeRUN)-thisBuff):(time2alz(timeRUN)+thisBuff));

	%detrending and zscoring
	for trialRUN = 1:1:size(dat1_t, 1)
		dat1_t(trialRUN, :) = zscore(detrend(dat1_t(trialRUN, :)));
		dat2_t(trialRUN, :) = zscore(detrend(dat2_t(trialRUN, :)));
	end %trialRUN = 1:1:size(dat1_t, 1)

	dat1_t = reshape(transpose(dat1_t), 1, nr * nl);
	dat2_t = reshape(transpose(dat2_t), 1, nr * nl);

	dat_t = [dat1_t ; dat2_t];

	[Ax Ex] = armorf(dat1_t, nr, nl, szOrder);
	[Ay Ey] = armorf(dat2_t, nr, nl, szOrder);
	[Axy E] = armorf(dat_t, nr, nl, szOrder);

	eyx = E(2, 2) - E(1, 2)^2 / E(1, 1);
	exy = E(1, 1) - E(2, 1)^2 / E(2, 2);
	N = size(E, 1);

	for freqRUN = 1:1:length(freq2use)
		H = eye(N);
		for m = 1:szOrder
			H = H + Axy(:,(m-1)*N+1:m*N)*exp(-1i*m*2*pi*freq2use(freqRUN)/SF);
		end	%m = 1:szOrder

		Hi = inv(H);
		S = H\E*Hi'/SF;

		GCMat(freqRUN, timeRUN, 1) = log( abs(S(2,2))/abs(S(2,2)-(Hi(2,1)*exy*conj(Hi(2,1)))/SF) );
		GCMat(freqRUN, timeRUN, 2) = log( abs(S(1,1))/abs(S(1,1)-(Hi(1,2)*eyx*conj(Hi(1,2)))/SF) );

		clear H Hi S;
	end	%freqRUN = 1:1:length(freq2use)

	clear dat1_t dat2_t dat_t Ax Ex Ay Ey Axy E eyx exy N;
end %timeRUN = 1:1:length(time2alz)

function varargout = armorf(x,Nr,Nl,p)
%ARMORF   AR parameter estimation via LWR method by Morf modified.
%   x is a matrix whose every row is one variable's time series
%   Nr is the number of realizations, Nl is the length of every realization
%   If the time series are stationary long, just let Nr=1, Nl=length(x)
%   p is the order of AR model
%
%   A = ARMORF(X,NR,NL,P) returns the polynomial coefficients A corresponding to 
%     the AR model estimate of matrix X using Morf's method.
%
%   [A,E] = ARMORF(...) returns the final prediction error E (the
%   covariance matrix of the white noise of the AR model).
%
%   [A,E,K] = ARMORF(...) returns the vector K of reflection 
%     coefficients (parcor coefficients).
%
%   Ref: M. Morf, etal, Recursive Multichannel Maximum Entropy Spectral Estimation,
%              IEEE trans. GeoSci. Elec., 1978, Vol.GE-16, No.2, pp85-94.
%        S. Haykin, Nonlinear Methods of Spectral Analysis, 2nd Ed.
%              Springer-Verlag, 1983, Chapter 2
%
%   finished on Aug.9, 2002 by Yonghong Chen

% Initialization
[L,N]=size(x);
R0=zeros(L,L);
R0f=R0;
R0b=R0;
pf=R0;
pb=R0;
pfb=R0;
ap(:,:,1)=R0;
bp(:,:,1)=R0;
En=R0;
for i=1:Nr
    En=En+x(:,(i-1)*Nl+1:i*Nl)*x(:,(i-1)*Nl+1:i*Nl)';
    ap(:,:,1)=ap(:,:,1)+x(:,(i-1)*Nl+2:i*Nl)*x(:,(i-1)*Nl+2:i*Nl)';        
    bp(:,:,1)=bp(:,:,1)+x(:,(i-1)*Nl+1:i*Nl-1)*x(:,(i-1)*Nl+1:i*Nl-1)';
end
ap(:,:,1) = inv((chol(ap(:,:,1)/Nr*(Nl-1)))');
bp(:,:,1) = inv((chol(bp(:,:,1)/Nr*(Nl-1)))');
for i=1:Nr
    efp = ap(:,:,1)*x(:,(i-1)*Nl+2:i*Nl);
    ebp = bp(:,:,1)*x(:,(i-1)*Nl+1:i*Nl-1);
    pf = pf + efp*efp';
    pb = pb + ebp*ebp';
    pfb = pfb + efp*ebp';
end
En = chol(En/N)'; % Covariance of the noise

% Initial output variables
coeff = [];%  Coefficient matrices of the AR model
kr=[];  % reflection coefficients

for m=1:p
   % Calculate the next order reflection (parcor) coefficient
   ck = inv((chol(pf))')*pfb*inv(chol(pb));
   kr=[kr,ck];
   % Update the forward and backward prediction errors
   ef = eye(L)- ck*ck';
   eb = eye(L)- ck'*ck;
     
   % Update the prediction error
   En = En*chol(ef)';
   E = (ef+eb)./2;   
   
   % Update the coefficients of the forward and backward prediction errors
   ap(:,:,m+1) = zeros(L);
   bp(:,:,m+1) = zeros(L);
   pf = zeros(L);
   pb = zeros(L);
   pfb = zeros(L);

   for i=1:m+1       
       a(:,:,i) = inv((chol(ef))')*(ap(:,:,i)-ck*bp(:,:,m+2-i));
       b(:,:,i) = inv((chol(eb))')*(bp(:,:,i)-ck'*ap(:,:,m+2-i));
   end
   for k=1:Nr
       efp = zeros(L,Nl-m-1);
       ebp = zeros(L,Nl-m-1);
       for i=1:m+1
           k1=m+2-i+(k-1)*Nl+1;
           k2=Nl-i+1+(k-1)*Nl;
           efp = efp+a(:,:,i)*x(:,k1:k2);
           ebp = ebp+b(:,:,m+2-i)*x(:,k1-1:k2-1);
       end
       pf = pf + efp*efp';
       pb = pb + ebp*ebp';
       pfb = pfb + efp*ebp';
   end
   ap = a;
   bp = b;
end
for j=1:p
    coeff = [coeff,inv(a(:,:,1))*a(:,:,j+1)];
end

varargout{1} = coeff;
if nargout >= 2
    varargout{2} = En*En';
end
if nargout >= 3
    varargout{3} = kr;
end
end
end	%jkComputeGranger(sig1, sig2, SF, szBuff, szSW, FreqRange)