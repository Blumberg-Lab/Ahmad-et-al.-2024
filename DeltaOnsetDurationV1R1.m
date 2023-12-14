%Sleep project by Blumberg lab
%Midha, Greta, & Mark
%
%Jangjin Kim, 2023-Aug-16

%initialization
clear all; close all; fclose all; clc;

%def basic paths
datROOT = ['\\lc-rs-store21.hpc.uiowa.edu\Blumberg_Lab_LSS\Jin\Data\'];
alzROOT = ['G:\Blumberg\PZProjectV3']; if ~exist(alzROOT) mkdir(alzROOT); end
	preprocROOT = [alzROOT '\preprocessV0']; if ~exist(preprocROOT) mkdir(preprocROOT); end
	grROOT = [alzROOT '\DeltaOnsetDurationV1R1']; if ~exist(grROOT) mkdir(grROOT); end

ageGROUP = {['P10'] ; ['P12']};
sleepStage = {['Active'] ; ['Quiet'] ; ['Wake']};

%basic params
SFclfp = 976.56;        %976.56Hz; since time interval is .001, I will regard this as 1000Hz
SFspk = 24.4 * 10^3;    %24.4kHz
SFResp = 30.518;        %30.518Hz

SFtarg = 1000;

%spectral param
roiFreq = [eps 30]; dcFreq = [eps .5]; grangerROIFreq = [eps 20];
deltaFreq= [.5 4];
thetaFreq = [4 7];

ratIDs2go = {['10D1'] ; ['12D2'] ; ['3J1'] ; ['6J3'] ; ['7B1'] ; ['7F1'] ; ['7M1'] ; ['7Q1']};
cohFlag = [1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1];
pzPLVFlag = [1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1];
m1PLVFlag = [1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1];
respFlag = [1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1];

txtOX = {['X'] ; ['O']}; txtDIR = {['PZ>M1'] ; ['M1>PZ']};
colLINEs = get(groot, 'DefaultAxesColorOrder'); colLINEs = [colLINEs ; 0 0 1]; colLINEs = [colLINEs ; fliplr(colLINEs)]; szTXT = 7.5;
if ~exist([grROOT '\deltaOnsetDurationTab.mat'])
	deltaOnsetDurationTab = [];

	for ageRUN = 1%2:1:size(ageGROUP, 1)
		for ratRUN = 1:1:size(ratIDs2go, 1)%size(ratIDs2go, 1):-1:1
		%for ratRUN = size(ratIDs2go, 1):-1:1
			thisRATID = ratIDs2go{ratRUN};
			disp(['Working on ' thisRATID]);

			preprocLOADROOT = [preprocROOT '\' ageGROUP{ageRUN} '\' thisRATID];
			load([preprocLOADROOT '\preproc.mat']);

			slTime = preproc.slTime;
			sResp = preproc.sResp;
			tSpace = preproc.tSpace;
			pzlfp = preproc.pzlfp;
			m1lfp = preproc.m1lfp;
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

			[pzDummy pzts2go] = jkParsePowerSection(pzDelta, lfpTSpace, SFtarg, 1, 3, 1); %1 med, 3 sec, + smoothing!
			[m1Dummy m1ts2go] = jkParsePowerSection(m1Delta, lfpTSpace, SFtarg, 1, 3, 1); %1 med, 3 sec, + smoothing!
            
            [pzDummy2 pzts2go_byPZ] = jkParsePowerSection(pzDelta(pzDelta1MedRaw == 1), lfpTSpace(pzDelta1MedRaw == 1), SFtarg, 1, 3, 1); %1 med, 3 sec, + smoothing!
            [m1Dummy2 m1ts2go_byPZ] = jkParsePowerSection(m1Delta(pzDelta1MedRaw == 1), lfpTSpace(pzDelta1MedRaw == 1), SFtarg, 1, 3, 1); %1 med, 3 sec, + smoothing!
            
            [pzDummy3 pzts2go_byM1] = jkParsePowerSection(pzDelta(m1Delta1MedRaw == 1), lfpTSpace(m1Delta1MedRaw == 1), SFtarg, 1, 3, 1); %1 med, 3 sec, + smoothing!
            [m1Dummy3 m1ts2go_byM1] = jkParsePowerSection(m1Delta(m1Delta1MedRaw == 1), lfpTSpace(m1Delta1MedRaw == 1), SFtarg, 1, 3, 1); %1 med, 3 sec, + smoothing!

			deltaOnsetDurationTab = [deltaOnsetDurationTab ; ... 
										[ageRUN ratRUN ... 
                                            size(pzts2go, 1) nanmean((pzts2go(:, 2) - pzts2go(:, 1)), 1) size(m1ts2go, 1) nanmean((m1ts2go(:, 2) - m1ts2go(:, 1)), 1) ... 
                                            size(pzts2go_byPZ, 1) nanmean((pzts2go_byPZ(:, 2) - pzts2go_byPZ(:, 1)), 1) size(m1ts2go_byPZ, 1) nanmean((m1ts2go_byPZ(:, 2) - m1ts2go_byPZ(:, 1)), 1) ... 
                                            size(pzts2go_byM1, 1) nanmean((pzts2go_byM1(:, 2) - pzts2go_byM1(:, 1)), 1) size(m1ts2go_byM1, 1) nanmean((m1ts2go_byM1(:, 2) - m1ts2go_byM1(:, 1)), 1)]];

			clear thisRATID preprocLOADROOT preproc slTime sResp tSpace pzlfp m1lfp lfpTSpace midhaAS midhaQS midhaWA midhaASraw midhaQSraw midhaWAraw islerQS islerAS islerQSraw islerASraw resp2SDraw resp2SD pzDelta m1Delta pzDelta1MedRaw pzDelta1Med m1Delta1MedRaw m1Delta1Med pzMUA pzSC pzRasters pzIFR grIndStruct;
            clear pzDummy pzts2go m1Dummy m1ts2go pzDummy2 pzts2go_byPZ m1Dummy2 m1ts2go_byPZ pzDummy3 pzts2go_byM1 m1Dummy3 m1ts2go_byM1;
		end %ratRUN = 1:1:size(ratIDs2go, 1)
	end %ageRUN = 2:1:size(ageGROUP, 1)
	save([grROOT '\deltaOnsetDurationTab.mat'], 'deltaOnsetDurationTab');
else
	load([grROOT '\deltaOnsetDurationTab.mat']);
end %~exist([grROOT '\grangerSpecTab.mat'])
sdtop

figure('Color', 'w', 'Position', [50 50 600 450]);
for measRUN = 1:1:2
	if measRUN == 1
		idx2go = 1:4;
	elseif measRUN == 2
		idx2go = (1:4) + 4;
	end %measRUN == 1

	subplot(1, 2, measRUN, 'FontSize', szTXT);
	boxplot(grangerMasterTab(cohFlag == 1, idx2go));
	set(gca, 'XTick', 1:1:4, 'XTickLabel', {['AS-PZ>M1'], ['AS-PZ<M1'], ['QS-PZ>M1'], ['QS-PZ<M1']})
end %measRUN = 1:1:2

nFIGROW = 1; nFIGCOL = 2;
picID = figure('Color', 'w', 'Position', [50 50 600 450]);
for slRUN = 1:1:2
	if slRUN == 1
		thisRatIdddx = 1:2:18; ln2gogo = 'b-'; ln2gogo2 = 'b:'; txtTitle2go = ['PZ > M1'];
	elseif slRUN == 2
		thisRatIdddx = 2:2:18; ln2gogo = 'r-'; ln2gogo2 = 'r:'; txtTitle2go = ['M1 > PZ'];
	end %slRUN == 1
	for dirRUN = 1:1:2
		if dirRUN == 1
			dirGOIdx = 1:20;
		elseif dirRUN == 2
			dirGOIdx = 21:40;
		end %dirRUN == 1
		thisRatDat2go = grangerSpecTab(thisRatIdddx, dirGOIdx);
			thisRatDat2go = thisRatDat2go(find(cohFlag == 1), :);

		subplot(nFIGROW, nFIGCOL, dirRUN, 'FontSize', szTXT);
		hold on;
			plot(freq2use, nanmean(thisRatDat2go, 1), ln2gogo, 'LineWidth', 2);
				plot(freq2use, nanmean(thisRatDat2go, 1) + nanstd(thisRatDat2go, [], 1) ./ sqrt(size(thisRatDat2go, 1)), ln2gogo2, 'LineWidth', .75);
				plot(freq2use, nanmean(thisRatDat2go, 1) - nanstd(thisRatDat2go, [], 1) ./ sqrt(size(thisRatDat2go, 1)), ln2gogo2, 'LineWidth', .75);
		hold off;
		set(gca, 'XLim', [.5 max(freq2use) + .5]); xlabel(['Frequency (Hz)']); ylabel(['G-connectivity']);
		if dirRUN == 1
			set(gca, 'YLim', [0 .6]);
		elseif dirRUN == 2
			set(gca, 'YLim', [0 .08]);
		end %dirRUN == 1
		title([txtTitle2go]);

		clear thisRatDat2go
	end %dirRUN = 1:1:2
end %slRUN = 1:1:2

stophere
baseFreq = [4 7];
deltaIdx2go = dsearchn(transpose(cohFreq2go), transpose(deltaFreq));
denomIdx2go = dsearchn(transpose(cohFreq2go), transpose(baseFreq));
wholeAS = cohMasterTab(cohFlag == 1, 1:50); wholeQS = cohMasterTab(cohFlag == 1, 51:100);

statMeas = []; statMeas2 = []; statMeas3 = [];
nFIGROW = 2; nFIGCOL = 4;
picID = figure('Color', 'w', 'Position', [50 50 600 450]);
for slRUN = 1:1:2
	if slRUN == 1
		thisWhole = wholeAS;
	elseif slRUN == 2
		thisWhole = wholeQS;
	end %slRUN == 1

	subplot(nFIGROW, nFIGCOL, [1 2] + (slRUN - 1) * 2, 'FontSize', szTXT);
	hold on;
		plot(cohFreq2go, nanmean(thisWhole, 1), '-k');
			plot(cohFreq2go, nanmean(thisWhole, 1) + nanstd(thisWhole, [], 1) ./ sqrt(nansum(cohFlag)), ':k');
			plot(cohFreq2go, nanmean(thisWhole, 1) - nanstd(thisWhole, [], 1) ./ sqrt(nansum(cohFlag)), ':k');
	hold off;
	set(gca, 'XTick', sort(unique([deltaFreq baseFreq])), 'YLim', [0 .8]); xlabel(['Frequency (Hz)']); if slRUN == 1 ylabel(['Coherence']); end
	title([sleepStage{slRUN} ' (n = ' num2str(nansum(cohFlag)) ')']);

	%compute measures
	thisMeas = nansum(thisWhole(:, deltaIdx2go(1):deltaIdx2go(2)), 2) ./ nansum(thisWhole(:, denomIdx2go(1):denomIdx2go(2)), 2); %ratio
	statMeas = [statMeas thisMeas];

	subplot(nFIGROW, nFIGCOL, [5 6], 'FontSize', szTXT);
	hold on;
		for indRUN = 1:1:size(thisWhole, 1)
			plot(slRUN, thisMeas(indRUN, 1), 'Color', ones(1, 3) .* .65, 'Marker', '.');
		end %indRUN = 1:1:size(thisWhole, 1)
		errorbar(slRUN, nanmean(thisMeas, 1), nanstd(thisMeas, [], 1) ./ sqrt(nansum(cohFlag)), 'ks');
	hold off;
	set(gca, 'XLim', [.5 2.5], 'XTick', [1 2], 'XTickLabel', {['AS'], ['QS']});
	ylabel(['detla / theta (' num2str(min(baseFreq)) '~' num2str(max(baseFreq)) 'Hz) coherence ratio']);

	%mean delta coh
	thisMeas2 = nanmean(thisWhole(:, deltaIdx2go(1):deltaIdx2go(2)), 2); %ratio
	statMeas2 = [statMeas2 thisMeas2];

	subplot(nFIGROW, nFIGCOL, [7], 'FontSize', szTXT);
	hold on;
		for indRUN = 1:1:size(thisWhole, 1)
			plot(slRUN, thisMeas2(indRUN, 1), 'Color', ones(1, 3) .* .65, 'Marker', '.');
		end %indRUN = 1:1:size(thisWhole, 1)
		errorbar(slRUN, nanmean(thisMeas2, 1), nanstd(thisMeas2, [], 1) ./ sqrt(nansum(cohFlag)), 'ks');
	hold off;
	set(gca, 'XLim', [.5 2.5], 'XTick', [1 2], 'XTickLabel', {['AS'], ['QS']});
	ylabel(['avg. detla coherence']);

	%coherence peak frequency
	thisMeas3 = [];
	for ratRUN = 1:1:size(thisWhole, 1)
		thisMeas3 = [thisMeas3 ; cohFreq2go(find(thisWhole(ratRUN, :) == nanmax(thisWhole(ratRUN, :), [], 2)))];
	end %ratRUN = 1:1:size(thisWhole, 1)
	statMeas3 = [statMeas3 thisMeas3];

	subplot(nFIGROW, nFIGCOL, [8], 'FontSize', szTXT);
	hold on;
		for indRUN = 1:1:size(thisWhole, 1)
			plot(slRUN, thisMeas3(indRUN, 1), 'Color', ones(1, 3) .* .65, 'Marker', '.');
		end %indRUN = 1:1:size(thisWhole, 1)
		errorbar(slRUN, nanmean(thisMeas3, 1), nanstd(thisMeas3, [], 1) ./ sqrt(nansum(cohFlag)), 'ks');
	hold off;
	set(gca, 'XLim', [.5 2.5], 'XTick', [1 2], 'XTickLabel', {['AS'], ['QS']});
	ylabel(['avg. coherence peak frequency (Hz)']);
end %slRUN = 1:1:2

[h p ci stats1] = ttest2(statMeas(:, 1), statMeas(:, 2))
[h2 p2 ci2 stats2] = ttest2(statMeas2(:, 1), statMeas2(:, 2))
[h3 p3 ci3 stats3] = ttest2(statMeas3(:, 1), statMeas3(:, 2))