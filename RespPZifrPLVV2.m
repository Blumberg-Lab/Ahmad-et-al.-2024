%Sleep project by Blumberg lab
%Midha, Greta, & Mark
%
%Jangjin Kim, 2023-Sep-21
%preliminary code: RespPhaseLockAlzV0R1.m in Version1

%initialization
clear all; close all; fclose all; clc;

%def basic paths
datROOT = ['\\lc-rs-store21.hpc.uiowa.edu\Blumberg_Lab_LSS\Jin\Data\'];
alzROOT = ['G:\Blumberg\PZProjectV3']; if ~exist(alzROOT) mkdir(alzROOT); end
	preprocROOT = [alzROOT '\preprocessV0']; if ~exist(preprocROOT) mkdir(preprocROOT); end
	resppzROOT = [alzROOT '\RespPZifrV2']; if ~exist(resppzROOT) mkdir(resppzROOT); end

ageGROUP = {['P10'] ; ['P12']};
sleepStage = {['Active'] ; ['Quiet'] ; ['Wake']};

%basic params
SFclfp = 976.56;        %976.56Hz; since time interval is .001, I will regard this as 1000Hz
SFspk = 24.4 * 10^3;    %24.4kHz
SFResp = 30.518;        %30.518Hz

SFtarg = 1000;

%spectral param
roiFreq = [eps 30]; dcFreq = [eps .5];
deltaFreq= [.5 4];
thetaFreq = [4 7];

ratIDs2go = {['10D1'] ; ['12D2'] ; ['3J1'] ; ['6J3'] ; ['7B1'] ; ['7F1'] ; ['7M1'] ; ['7Q1']};
cohFlag = [1 ; 1 ; 0 ; 1 ; 1 ; 1 ; 1 ; 1];
pzPLVFlag = [1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1];
m1PLVFlag = [1 ; 1 ; 0 ; 1 ; 1 ; 1 ; 1 ; 1];
respFlag = [1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1];

degSpace = 0:20:360; txtCRI22 = {['> 0 '] ; ['> .25'] ; ['> .5'] ; ['> .75'] ; ['> .9']};

txtOX = {['X'] ; ['O']};
colLINEs = get(groot, 'DefaultAxesColorOrder'); colLINEs = [colLINEs ; 0 0 1]; colLINEs = [colLINEs ; fliplr(colLINEs)]; szTXT = 7.5;
if ~exist([resppzROOT '\respPlvTab.mat'])
	respPlvTab = []; %50 x 2; first 50 for AS & secon 50 for QS

	for ageRUN = 1%2:1:size(ageGROUP, 1)
		for ratRUN = 1:1:size(ratIDs2go, 1)
			thisRATID = ratIDs2go{ratRUN};
			disp(['Working on ' thisRATID]);

			preprocLOADROOT = [preprocROOT '\' ageGROUP{ageRUN} '\' thisRATID];
			load([preprocLOADROOT '\preproc.mat']);

			slTime = preproc.slTime;
			sResp = preproc.sResp;
			tSpace = preproc.tSpace;
			pzlfp = preproc.pzlfp;
				pzlfplsf = downsample(pzlfp, round(length(pzlfp) / length(sResp)));
			m1lfp = preproc.m1lfp;
				m1lfplsf = downsample(m1lfp, round(length(pzlfp) / length(sResp)));
			lfpTSpace = preproc.lfpTSpace;
			midhaAS = preproc.midhaAS;
			midhaQS = preproc.midhaQS;
			midhaWA = preproc.midhaWA;
			midhaASraw = preproc.midhaASraw;
			midhaQSraw = preproc.midhaQSraw;
			midhaWAraw = preproc.midhaWAraw;
			islerQS = preproc.islerQS;
			islerAS = preproc.islerAS;
			islerQSraw = preproc.islerQSraw;
			islerASraw = preproc.islerASraw;
			resp2SDraw = preproc.resp2SDraw;
			resp2SD = preproc.resp2SD;
			pzDelta = preproc.pzDelta;
			m1Delta = preproc.m1Delta;
			pzDelta1MedRaw = preproc.pzDelta1MedRaw;
			pzDelta1Med = preproc.pzDelta1Med;
			m1Delta1MedRaw = preproc.m1Delta1MedRaw;
			m1Delta1Med = preproc.m1Delta1Med;
			pzMUA = preproc.pzMUA;
			pzSC = preproc.pzSC;
			pzRasters = preproc.pzRasters;
			pzIFR = preproc.pzIFR;
				smpzIFR = smooth(pzIFR);
				envpzIFR = envelope(smpzIFR, 500, 'analytic');
				smpzIFRAbvThr = nan(length(smpzIFR), 1);
				smpzIFRAbvThr(find(smpzIFR > (median(envpzIFR) * 1))) = smpzIFR(find(smpzIFR > (median(envpzIFR) * 1)));

			%get phase info of respiratory signals
			phaseResp = angle(hilbert(sResp));
			dspzIFR = downsample(smpzIFRAbvThr, round(length(smpzIFRAbvThr) / length(phaseResp)));
		
			pzifrphase2put = []; thisYLim = [];
			%QS only
			nFIGROW = 2; nFIGCOL = 5;
			picID = figure('Color', 'w', 'Position', [25 25 1250 600]);
			for criRUN = 1:1:5 	%0, .25, .5, .75, .9
				dspzIFRidx = zeros(1, length(phaseResp));
				if criRUN == 1
					dspzIFRidx(find(dspzIFR > 0)) = 1; %dspzIFRidx(find(dspzIFR > nanmean(dspzIFR))) = 1;
				elseif criRUN == 2
					dspzIFRidx(find(dspzIFR > .25)) = 1; %dspzIFRidx(find(dspzIFR > nanmean(dspzIFR))) = 1;
				elseif criRUN == 3
					dspzIFRidx(find(dspzIFR > .5)) = 1; %dspzIFRidx(find(dspzIFR > nanmean(dspzIFR))) = 1;
				elseif criRUN == 4
					dspzIFRidx(find(dspzIFR > .75)) = 1; %dspzIFRidx(find(dspzIFR > nanmean(dspzIFR))) = 1;
				elseif criRUN == 5
					dspzIFRidx(find(dspzIFR > .9)) = 1; %dspzIFRidx(find(dspzIFR > nanmean(dspzIFR))) = 1;
				end %criRUN == 1

				%stateDefIFR = [resp2SD ; midhaQS ; islerQS ; pzDelta1Med ; m1Delta1Med ; dspzIFRidx];
                stateDefIFR = [resp2SD ; midhaQS ; islerQS ; pzDelta1Med ; dspzIFRidx];
				dat2go = find(nansum(stateDefIFR, 1) == size(stateDefIFR, 1));

				subplot(nFIGROW, nFIGCOL, criRUN, 'FontSize', szTXT);
				hold on;
					bar(degSpace, histc(mod(rad2deg(phaseResp(dat2go)) + 720, 360), degSpace));
				hold off;
				set(gca, 'XLim', [-.5 360.5], 'XTick', [0 180 360]);
				title({[ratIDs2go{ratRUN} '-' sleepStage{2}] ; [txtCRI22{criRUN}]});

				subplot(nFIGROW, nFIGCOL, criRUN + nFIGCOL, 'FontSize', szTXT);
				polarhistogram(phaseResp(dat2go), length(degSpace));
				thetaticks([0 180 360]);

				if length(dat2go) == 0
					thisHist2go = nan(1, length(degSpace));
				else
					thisHist2go = histc(mod(rad2deg(phaseResp(dat2go)) + 720, 360), degSpace) ./ nansum(histc(mod(rad2deg(phaseResp(dat2go)) + 720, 360), degSpace));
					if size(thisHist2go, 1) > size(thisHist2go, 2)
						thisHist2go = transpose(thisHist2go);
					end %size(thisHist2go, 1) > size(thisHist2go, 2)
				end %length(dat2go) == 0

				pzifrphase2put = [pzifrphase2put thisHist2go(1:18)];

				clear stateDefIFR dat2go thisHist2go;
			end %criRUN = 1:1:4 	%0, .5, .75, .9

			saveas(picID, [resppzROOT '\' ratIDs2go{ratRUN} '-RespPZIfr.bmp']); close all;
			respPlvTab = [respPlvTab ; pzifrphase2put];

			clear thisRATID preprocLOADROOT preproc slTime sResp tSpace pzlfp m1lfp lfpTSpace midhaAS midhaQS midhaWA midhaASraw midhaQSraw midhaWAraw islerQS islerAS islerQSraw islerASraw resp2SDraw resp2SD pzDelta m1Delta pzDelta1MedRaw pzDelta1Med m1Delta1MedRaw m1Delta1Med pzMUA pzSC pzRasters pzIFR;
		end %ratRUN = 1:1:size(ratIDs2go, 1)
	end %ageRUN = 2:1:size(ageGROUP, 1)
	save([resppzROOT '\respPlvTab.mat'], 'respPlvTab');
else
	load([resppzROOT '\respPlvTab.mat']);
end %~exist([resppzROOT '\respPlvTab.mat'])

qsRPZIFR_0_1 = 1; qsRPZIFR_0_18 = qsRPZIFR_0_1 + 17;
qsRPZIFR_quart_1 = qsRPZIFR_0_18 + 1; qsRPZIFR_quart_18 = qsRPZIFR_quart_1 + 17;
qsRPZIFR_half_1 = qsRPZIFR_quart_18 + 1; qsRPZIFR_half_18 = qsRPZIFR_half_1 + 17; 
qsRPZIFR_3quart_1 = qsRPZIFR_half_18 + 1; qsRPZIFR_3quart_18 = qsRPZIFR_3quart_1 + 17;
qsRPZIFR_high_1 = qsRPZIFR_3quart_18 + 1; qsRPZIFR_high_18 = qsRPZIFR_high_1 + 17;

nFIGROW = 1; nFIGCOL = 2; xaxisVals = 10:20:350; tab4rawVal = [];
criRUN = 1; 

picID = figure('Color', 'w', 'Position', [25 25 1250 400]);
if criRUN == 1
	sttIdxx = qsRPZIFR_0_1;
elseif criRUN == 2
	sttIdxx = qsRPZIFR_quart_1;
elseif criRUN == 3
	sttIdxx = qsRPZIFR_half_1;
end 
eddIdxx = sttIdxx + 17;
nRats2use = nansum(pzPLVFlag);

smFac = 3;
thisAvg = nanmean(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), 1);

subplot(nFIGROW, nFIGCOL, 1, 'FontSize', szTXT);
hold on;
	%errorbar(10:20:350, nanmean(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), 1), nanstd(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), [], 1) ./ sqrt(nRats2use), '-bs');
	% plot(xaxisVals, nanmean(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), 1), '-b', 'LineWidth', 2);
	% 	plot(xaxisVals, nanmean(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), 1) + nanstd(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), [], 1) ./ sqrt(nRats2use), ':b', 'LineWidth', 1);
	% 	plot(xaxisVals, nanmean(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), 1) - nanstd(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), [], 1) ./ sqrt(nRats2use), ':b', 'LineWidth', 1);
	bar([10:20:710], [smooth(thisAvg, smFac) ; smooth(thisAvg, smFac)]);
hold off;
set(gca, 'XLim', [5 715], 'XTick', [170 330], 'XTickLabel', {['T'], ['P']}); title({['Pop resp-PZ ifr plv (N = ' num2str(nRats2use) ')'] ; [txtCRI22{criRUN}]});
if criRUN <= 3 set(gca, 'YLim', [0 .15]); else set(gca, 'YLim', [0 .3]); end

tab4rawVal = [tab4rawVal respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx)];

thisAvg = smooth(thisAvg, smFac);
subplot(nFIGROW, nFIGCOL, 2, 'FontSize', szTXT);
hold on;
	for drawRUN = 1:1:length(xaxisVals)
		thisTheta = linspace((drawRUN - 1) * 20, (drawRUN - .1) * 20, 21);
		patch([0 thisAvg(drawRUN) * cosd(thisTheta)], [0 thisAvg(drawRUN) * sind(thisTheta)], 'b'); alpha(.75);

		if xaxisVals(drawRUN) == 10 | xaxisVals(drawRUN) == 170 | xaxisVals(drawRUN) == 330
			rtext = thisAvg(drawRUN) * 1.15;
    		text(rtext * cosd((drawRUN - .5) * 20), rtext * sind((drawRUN - .5) * 20), num2str(xaxisVals(drawRUN)));
    	end %xaxisVals(drawRUN) == 10 | xaxisVals(drawRUN) == 170 | xaxisVals(drawRUN) == 330
	end %drawRUN = 1:1:length(xaxisVals)
	%thetaticks([0 175 330]);
hold off;
axis equal; axis on;

% nFIGROW = 2; nFIGCOL = 4; xaxisVals = 10:20:350; tab4rawVal = [];
% picID = figure('Color', 'w', 'Position', [25 25 1250 400]);
% for criRUN = 1:1:nFIGCOL
% 	if criRUN == 1
% 		sttIdxx = qsRPZIFR_0_1;
% 	elseif criRUN == 2
% 		sttIdxx = qsRPZIFR_quart_1;
% 	elseif criRUN == 3
% 		sttIdxx = qsRPZIFR_half_1;
% 	elseif criRUN == 4
% 		sttIdxx = qsRPZIFR_3quart_1;
% 	elseif criRUN == 5
% 		sttIdxx = qsRPZIFR_high_1;
% 	end %criRUN == 1
% 	eddIdxx = sttIdxx + 17;
% 	nRats2use = nansum(pzPLVFlag);

% 	subplot(nFIGROW, nFIGCOL, criRUN, 'FontSize', szTXT);
% 	hold on;
% 		%errorbar(10:20:350, nanmean(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), 1), nanstd(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), [], 1) ./ sqrt(nRats2use), '-bs');
% 		plot(xaxisVals, nanmean(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), 1), '-b', 'LineWidth', 2);
% 			plot(xaxisVals, nanmean(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), 1) + nanstd(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), [], 1) ./ sqrt(nRats2use), ':b', 'LineWidth', 1);
% 			plot(xaxisVals, nanmean(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), 1) - nanstd(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), [], 1) ./ sqrt(nRats2use), ':b', 'LineWidth', 1);
% 	hold off;
% 	set(gca, 'XLim', [5 355], 'XTick', [175 330], 'XTickLabel', {['T'], ['P']}); title({['Pop resp-PZ ifr plv (N = ' num2str(nRats2use) ')'] ; [txtCRI22{criRUN}]});
% 	if criRUN <= 3 set(gca, 'YLim', [0 .15]); else set(gca, 'YLim', [0 .3]); end

% 	tab4rawVal = [tab4rawVal respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx)];

% 	thisAvg = nanmean(respPlvTab(find(pzPLVFlag == 1), sttIdxx:eddIdxx), 1);
% 	subplot(nFIGROW, nFIGCOL, criRUN + nFIGCOL, 'FontSize', szTXT);
% 	hold on;
% 		for drawRUN = 1:1:length(xaxisVals)
% 			thisTheta = linspace((drawRUN - 1) * 20, (drawRUN - .1) * 20, 21);
% 			patch([0 thisAvg(drawRUN) * cosd(thisTheta)], [0 thisAvg(drawRUN) * sind(thisTheta)], 'b'); alpha(.75);

% 			if xaxisVals(drawRUN) == 10 | xaxisVals(drawRUN) == 170 | xaxisVals(drawRUN) == 330
% 				rtext = thisAvg(drawRUN) * 1.15;
% 	    		text(rtext * cosd((drawRUN - .5) * 20), rtext * sind((drawRUN - .5) * 20), num2str(xaxisVals(drawRUN)));
% 	    	end %xaxisVals(drawRUN) == 10 | xaxisVals(drawRUN) == 170 | xaxisVals(drawRUN) == 330
% 		end %drawRUN = 1:1:length(xaxisVals)
% 		%thetaticks([0 175 330]);
% 	hold off;
% 	axis equal; axis on;
% end %criRUN = 1:1:5