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
	respROOT = [alzROOT '\RespAlignedV1R2']; if ~exist(respROOT) mkdir(respROOT); end

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
cohFlag = [1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1];
pzPLVFlag = [1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1];
m1PLVFlag = [1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1];
respFlag = [1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1];

txtOX = {['X'] ; ['O']}; txtLFP = {['PZ'] ; ['M1']};
colLINEs = get(groot, 'DefaultAxesColorOrder'); colLINEs = [colLINEs ; 0 0 1]; colLINEs = [colLINEs ; fliplr(colLINEs)]; szTXT = 7.5;
if ~exist([respROOT '\RespAligned.mat'])
	RespAligned = []; 

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
				pzDeltalsf = downsample(pzDelta, round(length(pzDelta) / length(sResp)));
				pzDeltaM = nanmean(abs(pzDeltalsf) .^ 2); pzDeltaSD = nanstd(abs(pzDeltalsf) .^ 2);
			m1Delta = preproc.m1Delta;
				m1Deltalsf = downsample(m1Delta, round(length(m1Delta) / length(sResp)));
				m1DeltaM = nanmean(abs(m1Deltalsf) .^ 2); m1DeltaSD = nanstd(abs(m1Deltalsf) .^ 2);
			pzDelta1MedRaw = preproc.pzDelta1MedRaw;
			pzDelta1Med = preproc.pzDelta1Med;
			m1Delta1MedRaw = preproc.m1Delta1MedRaw;
			m1Delta1Med = preproc.m1Delta1Med;
			pzMUA = preproc.pzMUA;
			pzSC = preproc.pzSC;
			pzRasters = preproc.pzRasters;
			pzIFR = preproc.pzIFR;

			nFIGROW = 2; nFIGCOL = 2;
			picID = figure('Color', 'w', 'Position', [50 50 650 600], 'Visible', 0);
			for slRUN = 1:1:2 	%AS & QS
				if slRUN == 1
					stateDef = [resp2SD ; midhaAS ; islerAS];
				elseif slRUN == 2
					%stateDef = [resp2SD ; midhaQS ; islerQS ; pzDelta1Med ; m1Delta1Med];
					stateDef = [resp2SD ; midhaQS ; islerQS ; pzDelta1Med];
				end %slRUN == 1
				dat2go = find(nansum(stateDef, 1) == size(stateDef, 1));

				thisSResp = sResp(dat2go); thisPZ = pzlfplsf(dat2go); thisM1 = m1lfplsf(dat2go);
				thisPZdel = pzDeltalsf(dat2go); thisM1del = m1Deltalsf(dat2go);
					thisPZdelPow = (abs(thisPZdel) .^ 2 - pzDeltaM) ./ pzDeltaSD;
					thisM1delPow = (abs(thisM1del) .^ 2 - m1DeltaM) ./ m1DeltaSD;

				thisPhaseSResp = mod(rad2deg(angle(hilbert(thisSResp))) + 360, 360);
				thisPZM1phasePow = nan(2, 360);
				for phaseRUN = 1:1:length(thisPhaseSResp)
					thisPZM1phasePow(1, floor(thisPhaseSResp(phaseRUN)) + 1) = nanmean([thisPZM1phasePow(1, floor(thisPhaseSResp(phaseRUN)) + 1) ... 
																						thisPZdelPow(phaseRUN)], 2);
					thisPZM1phasePow(2, floor(thisPhaseSResp(phaseRUN)) + 1) = nanmean([thisPZM1phasePow(2, floor(thisPhaseSResp(phaseRUN)) + 1) ... 
																						thisM1delPow(phaseRUN)], 2);
				end %phaseRUN = 1:1:length(thisPhaseSResp)

				for lfpRUN = 1:1:2
					subplot(nFIGROW, nFIGCOL, lfpRUN + (slRUN - 1) * nFIGCOL, 'FontSize', szTXT);
					hold on;
						plot(thisPZM1phasePow(lfpRUN, :));
						plot(smooth(thisPZM1phasePow(lfpRUN, :), 20));
					hold off;
					set(gca, 'XLim', [.5 360.6], 'XTick', [180]);
					xlabel(['Respiratory phase']); ylabel(['Avg. Z-Delta power']); title({[sleepStage{slRUN}] ; ['RespPhase-' txtLFP{lfpRUN} ' delta power']});

					RespAligned{ratRUN, lfpRUN + (slRUN - 1) * 2} = thisPZM1phasePow(lfpRUN, :);
				end %lfpRUN = 1:1:2
			end %slRUN = 1:1:2
			saveas(picID, [respROOT '\' ratIDs2go{ratRUN} '.bmp']); close all;

			clear thisRATID preprocLOADROOT preproc slTime sResp tSpace pzlfp m1lfp lfpTSpace midhaAS midhaQS midhaWA midhaASraw midhaQSraw midhaWAraw islerQS islerAS islerQSraw islerASraw resp2SDraw resp2SD pzDelta m1Delta pzDelta1MedRaw pzDelta1Med m1Delta1MedRaw m1Delta1Med pzMUA pzSC pzRasters pzIFR;
		end %ratRUN = 1:1:size(ratIDs2go, 1)
	end %ageRUN = 2:1:size(ageGROUP, 1)
	save([respROOT '\RespAligned.mat'], 'RespAligned');
else
	load([respROOT '\RespAligned.mat']);
end %~exist([respROOT '\RespAligned.mat'])

nFIGROW = 2; nFIGCOL = 2; xls2go = []; sm2use = 20;
picID = figure('Color', 'w', 'Position', [50 50 650 600]);
slRUN = 2;
for lfpRUN = 1:1:2
	if lfpRUN == 1
		thisFlg2go = pzPLVFlag;
	else
		thisFlg2go = m1PLVFlag;
	end %lfpRUN == 1
	ratTab = [];
	for ratRUN = 1:1:size(ratIDs2go, 1)
		ratTab = [ratTab ; ... 
					RespAligned{ratRUN, lfpRUN + (slRUN - 1) * 2}];
	end %ratRUN = 1:1:size(ratIDs2go, 1)
	xls2go{lfpRUN} = ratTab;
		ratTab(:, 360) = nanmean([ratTab(:, 360) ratTab(:, 1)], 2);

	subplot(nFIGROW, nFIGCOL, lfpRUN, 'FontSize', szTXT);
	hold on;
		plot(1:720, repmat(smooth(nanmean(ratTab(find(thisFlg2go == 1), :), 1), sm2use), 2, 1), 'Color', [1 0 0], 'LineStyle', '-', 'LineWidth', 1);
			plot(1:720, repmat(smooth(nanmean(ratTab(find(thisFlg2go == 1), :), 1) + nanstd(ratTab(find(thisFlg2go == 1), :), [], 1) ./ sqrt(nansum(thisFlg2go, 1)), sm2use), 2, 1), 'Color', [1 0 0], 'LineStyle', ':', 'LineWidth', .5);
			plot(1:720, repmat(smooth(nanmean(ratTab(find(thisFlg2go == 1), :), 1) - nanstd(ratTab(find(thisFlg2go == 1), :), [], 1) ./ sqrt(nansum(thisFlg2go, 1)), sm2use), 2, 1), 'Color', [1 0 0], 'LineStyle', ':', 'LineWidth', .5);				
	hold off;
	set(gca, 'XLim', [.5 360.5], 'XTick', [180 360 540], 'XTickLabel', {['pi'], ['2pi'], ['3pi']});
		set(gca, 'YLim', [-.5 2])
	if slRUN == 1 ylabel(['Avg. Z-delta power']); end
	title({[sleepStage{slRUN} ' (n = ' num2str(nansum(thisFlg2go, 1)) ')'] ; ['RespPhase-' txtLFP{lfpRUN} ' delta power']});

	subplot(nFIGROW, nFIGCOL, lfpRUN + nFIGCOL, 'FontSize', szTXT);
	%polarscatter(deg2rad(1:360), smooth(nanmean(ratTab(find(thisFlg2go == 1), :), 1), sm2use), 's')
	polarplot(deg2rad([1:360 1]), smooth(nanmean([ratTab(find(thisFlg2go == 1), :) ratTab(find(thisFlg2go == 1), 1)] , 1), sm2use), 'r-');
	set(gca, 'thetatick', [0 180], 'thetatickLabel', {['0, 2pi'], ['pi']});
		set(gca, 'RLim', [0 1.25], 'RTick', [.5 1])
end %lfpRUN = 1:1:2
jfjfas

nFIGROW = 2; nFIGCOL = 2; xls2go = []; sm2use = 20;
picID = figure('Color', 'w', 'Position', [50 50 650 600]);
slRUN = 2;
for lfpRUN = 1:1:2
	if lfpRUN == 1
		thisFlg2go = pzPLVFlag;
	else
		thisFlg2go = m1PLVFlag;
	end %lfpRUN == 1
	ratTab = [];
	for ratRUN = 1:1:size(ratIDs2go, 1)
		ratTab = [ratTab ; ... 
					RespAligned{ratRUN, lfpRUN + (slRUN - 1) * 2}];
	end %ratRUN = 1:1:size(ratIDs2go, 1)
	xls2go{lfpRUN} = ratTab;

	subplot(nFIGROW, nFIGCOL, lfpRUN, 'FontSize', szTXT);
	hold on;
		plot(1:720, repmat(smooth(nanmean(ratTab(find(thisFlg2go == 1), :), 1), sm2use), 2, 1), 'Color', [1 0 0], 'LineStyle', '-', 'LineWidth', 1);
			plot(1:720, repmat(smooth(nanmean(ratTab(find(thisFlg2go == 1), :), 1) + nanstd(ratTab(find(thisFlg2go == 1), :), [], 1) ./ sqrt(nansum(thisFlg2go, 1)), sm2use), 2, 1), 'Color', [1 0 0], 'LineStyle', ':', 'LineWidth', .5);
			plot(1:720, repmat(smooth(nanmean(ratTab(find(thisFlg2go == 1), :), 1) - nanstd(ratTab(find(thisFlg2go == 1), :), [], 1) ./ sqrt(nansum(thisFlg2go, 1)), sm2use), 2, 1), 'Color', [1 0 0], 'LineStyle', ':', 'LineWidth', .5);				
	hold off;
	set(gca, 'XLim', [.5 720.5], 'XTick', [180 360 540], 'XTickLabel', {['pi/2'], ['pi'], ['pi/2']});
	if slRUN == 1 ylabel(['Avg. Z-delta power']); end
	title({[sleepStage{slRUN} ' (n = ' num2str(nansum(thisFlg2go, 1)) ')'] ; ['RespPhase-' txtLFP{lfpRUN} ' delta power']});

	subplot(nFIGROW, nFIGCOL, lfpRUN + nFIGCOL, 'FontSize', szTXT);
	%polarscatter(deg2rad(1:360), smooth(nanmean(ratTab(find(thisFlg2go == 1), :), 1), sm2use), 's')
	polarplot(deg2rad([1:360 1]), smooth(nanmean([ratTab(find(thisFlg2go == 1), :) ratTab(find(thisFlg2go == 1), 1)] , 1), sm2use), 'r-');
	set(gca, 'thetatick', [0 180], 'thetatickLabel', {['0, pi'], ['pi/2']})
end %lfpRUN = 1:1:2
fasfas
nFIGROW = 2; nFIGCOL = 2; xls2go = []; sm2use = 20;
picID = figure('Color', 'w', 'Position', [50 50 650 600]);
for slRUN = 1:1:2
	for lfpRUN = 1:1:2
		if lfpRUN == 1
			thisFlg2go = pzPLVFlag;
		else
			thisFlg2go = m1PLVFlag;
		end %lfpRUN == 1
		ratTab = [];
		for ratRUN = 1:1:size(ratIDs2go, 1)
			ratTab = [ratTab ; ... 
						RespAligned{ratRUN, lfpRUN + (slRUN - 1) * 2}];
		end %ratRUN = 1:1:size(ratIDs2go, 1)
		xls2go{lfpRUN + (slRUN - 1) * 2} = ratTab;

		subplot(nFIGROW, nFIGCOL, slRUN + (lfpRUN - 1) * nFIGCOL, 'FontSize', szTXT);
		hold on;
			plot(1:360, nanmean(ratTab(find(thisFlg2go == 1), :), 1), 'Color', ones(1, 3) .* .65, 'LineStyle', '-', 'LineWidth', 1);
				plot(1:360, nanmean(ratTab(find(thisFlg2go == 1), :), 1) + nanstd(ratTab(find(thisFlg2go == 1), :), [], 1) ./ sqrt(nansum(thisFlg2go, 1)), 'Color', ones(1, 3) .* .65, 'LineStyle', ':', 'LineWidth', .5);
				plot(1:360, nanmean(ratTab(find(thisFlg2go == 1), :), 1) - nanstd(ratTab(find(thisFlg2go == 1), :), [], 1) ./ sqrt(nansum(thisFlg2go, 1)), 'Color', ones(1, 3) .* .65, 'LineStyle', ':', 'LineWidth', .5);
			plot(1:360, smooth(nanmean(ratTab(find(thisFlg2go == 1), :), 1), sm2use), 'Color', [1 0 0], 'LineStyle', '-', 'LineWidth', 1);
				plot(1:360, smooth(nanmean(ratTab(find(thisFlg2go == 1), :), 1) + nanstd(ratTab(find(thisFlg2go == 1), :), [], 1) ./ sqrt(nansum(thisFlg2go, 1)), sm2use), 'Color', [1 0 0], 'LineStyle', ':', 'LineWidth', .5);
				plot(1:360, smooth(nanmean(ratTab(find(thisFlg2go == 1), :), 1) - nanstd(ratTab(find(thisFlg2go == 1), :), [], 1) ./ sqrt(nansum(thisFlg2go, 1)), sm2use), 'Color', [1 0 0], 'LineStyle', ':', 'LineWidth', .5);				
		hold off;
		set(gca, 'XLim', [.5 360.6], 'XTick', [180]);
		if slRUN == 1 ylabel(['Avg. Z-delta power']); end
		title({[sleepStage{slRUN} ' (n = ' num2str(nansum(thisFlg2go, 1)) ')'] ; ['RespPhase-' txtLFP{lfpRUN} ' delta power']});
	end %lfpRUN = 1:1:2
end %slRUN = 1:1:2	

fafsa
ratTab = [];
nFIGROW = 2; nFIGCOL = 2; txtLFPs = {['PZ'] ; ['M1']};
for ratRUN = 1:1:size(ratIDs2go, 1)
	val2put = [];
	picID = figure('Color', 'w', 'Position', [45 45 600 600]); val2put = [];
	for slRUN = 1:1:2
		for lfpRUN = 1:1:2
			thisDat2vis = RespAligned{ratRUN, lfpRUN + (slRUN - 1) * 2};
				val2put = [val2put nanmean(thisDat2vis, 1)];

			subplot(nFIGROW, nFIGCOL, lfpRUN + (slRUN - 1) * nFIGCOL, 'FontSize', szTXT);
			hold on;
				%boxplot(thisDat2vis);
				errorbar(1:1:2, nanmean(thisDat2vis, 1), ... 
							nanstd(thisDat2vis, [], 1) ./ sqrt(size(thisDat2vis, 1)));
				val2put = [val2put nanmean(thisDat2vis, 1)];
			hold off;
			set(gca, 'XLim', [.5 2.5], 'XTick', 1:1:2, 'XTickLabel', {['inh.'], ['exh.']});
			ylabel(['z-delta power']);
			title({[ratIDs2go{ratRUN}] ; [sleepStage{slRUN} '-' txtLFPs{lfpRUN} ' (n = ' num2str(size(thisDat2vis, 1)) ')']});

			clear thisDat2vis;
		end %lfpRUN = 1:1:2
	end %slRUN = 1:1:2
	ratTab = [ratTab ; val2put];
end %ratRUN = 1:1:size(ratIDs2go, 1)

picID = figure('Color', 'w', 'Position', [45 45 600 600]); thisMeas = [];
for slRUN = 1:1:2
	for lfpRUN = 1:1:2
		if lfpRUN == 1
			thisFlag = pzPLVFlag;
		elseif lfpRUN == 2
			thisFlag = m1PLVFlag;
		end %lfpRUN == 1

		thisDat2vis = ratTab(thisFlag == 1, ([1 2] + (lfpRUN - 1) * 2 + (slRUN - 1) * 4));
		thisMeas{lfpRUN + (slRUN - 1) * 2} = thisDat2vis;

		subplot(nFIGROW, nFIGCOL, lfpRUN + (slRUN - 1) * nFIGCOL, 'FontSize', szTXT);
		hold on;
			boxplot(thisDat2vis);
		hold off;
		set(gca, 'XLim', [.5 2.5], 'XTick', 1:1:2, 'XTickLabel', {['inh.'], ['exh.']});
		ylabel(['z-delta power']);
		title({['Pop (n = ' num2str(size(thisDat2vis, 1)) ')'] ; [sleepStage{slRUN} '-' txtLFPs{lfpRUN}]});

		clear thisDat2vis;
	end %lfpRUN = 1:1:2
end %slRUN = 1:1:2

[h1 p1 ci1 stats1] = ttest2(thisMeas{1}(:, 1), thisMeas{1}(:, 2))
[h2 p2 ci2 stats2] = ttest2(thisMeas{2}(:, 1), thisMeas{2}(:, 2))
[h3 p3 ci3 stats3] = ttest2(thisMeas{3}(:, 1), thisMeas{3}(:, 2))
[h4 p4 ci4 stats4] = ttest2(thisMeas{4}(:, 1), thisMeas{4}(:, 2))

[h1 p1 ci1 stats1] = ttest(thisMeas{1}(:, 1), thisMeas{1}(:, 2))
[h2 p2 ci2 stats2] = ttest(thisMeas{2}(:, 1), thisMeas{2}(:, 2))
[h3 p3 ci3 stats3] = ttest(thisMeas{3}(:, 1), thisMeas{3}(:, 2))
[h4 p4 ci4 stats4] = ttest(thisMeas{4}(:, 1), thisMeas{4}(:, 2))