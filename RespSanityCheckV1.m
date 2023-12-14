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
	respROOT = [alzROOT '\RespSanityV1']; if ~exist(respROOT) mkdir(respROOT); end

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

txtOX = {['X'] ; ['O']};
colLINEs = get(groot, 'DefaultAxesColorOrder'); colLINEs = [colLINEs ; 0 0 1]; colLINEs = [colLINEs ; fliplr(colLINEs)]; szTXT = 7.5;
if ~exist([respROOT '\respTab.mat'])
	respTab = []; %50 x 2; first 50 for AS & secon 50 for QS

	for ageRUN = 1:1:1%size(ageGROUP, 1)
		for ratRUN = 1:1:size(ratIDs2go, 1)
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
			pzMUA = preproc.pzMUA;
			pzSC = preproc.pzSC;
			pzRasters = preproc.pzRasters;
			pzIFR = preproc.pzIFR;

			%get regularity of the respiratory signal
			[thisPks thisPksLocs] = findpeaks(sResp, 'MinPeakDistance', .25, 'MinPeakProminence', .00001); %thisPksLocs in s resolution
			RespPksIdx = zeros(1, length(sResp)); RespPksIdx(thisPksLocs) = 1;

			sRespWAripped = sResp; sRespWAripped(midhaWA == 1) = nan;
			[islerPks islerPksLocs] = findpeaks(sRespWAripped, tSpace, 'MinPeakDistance', .25, 'MinPeakProminence', .00001); %thisPksLocs in s resolution
			b2bInterval = diff(islerPksLocs);
			thisIBR = ones(length(b2bInterval), 1) ./ b2bInterval; %thisIBR(thisIBR < .5) = nan;
			thisIBR = [nan ; thisIBR]; varMat = nan(1, length(tSpace));
			for tRUN = 1:1:length(tSpace)
				tIdx2use = dsearchn(islerPksLocs, tSpace(tRUN));
				if abs(islerPksLocs(tIdx2use) - tSpace(tRUN)) < 1
					minTIdx = dsearchn(islerPksLocs, (tSpace(tRUN) - 1));
					maxTIdx = dsearchn(islerPksLocs, (tSpace(tRUN) + 1));

					period2use = [thisIBR(minTIdx:maxTIdx)];
					varMat(tRUN) = nanvar(period2use);
				end %abs(islerPksLocs(tIdx2use) - tSpace(tRUN)) < 1
				clear tIdx2use minTIdx maxTIdx period2use;
			end %tRUN = 1:1:length(tSpace)
			zVarMat = (varMat - nanmean(varMat, 2)) ./ nanstd(varMat);

			resp2put = [];
			nFIGROW = 1; nFIGCOL = 4;
			picID = figure('Color', 'w', 'Position', [25 25 890 410]);
			for slRUN = 1:1:2 	%AS & QS
				if slRUN == 1
					stateDef = [resp2SD ; midhaAS ; islerAS];
				elseif slRUN == 2
					%stateDef = [resp2SD ; midhaQS ; islerQS ; pzDelta1Med ; m1Delta1Med];
					stateDef = [resp2SD ; midhaQS ; islerQS ; pzDelta1Med];
				end %slRUN == 1
				dat2go = find(nansum(stateDef, 1) == size(stateDef, 1));

				%sleep duration
				subplot(nFIGROW, nFIGCOL, 1, 'FontSize', szTXT);
				hold on;
					plot(slRUN, length(dat2go) * SFResp / 1000 / 60, 'ks');
				hold off;
				if slRUN == 2 set(gca, 'XLim', [.5 2.5], 'XTick', 1:1:2, 'XTickLabel', {['AS'], ['QS']}); title({[ratIDs2go{ratRUN}] ; ['Resp duration (m)']}); end

				%peak interval
				subplot(nFIGROW, nFIGCOL, 2, 'FontSize', szTXT);
				hold on;
					errorbar(slRUN, nanmean(diff(find(RespPksIdx(dat2go) == 1))) .* SFResp, ... 
								nanstd(diff(find(RespPksIdx(dat2go) == 1)), [], 2) .* SFResp, 'ks');
				hold off;
				if slRUN == 2 set(gca, 'XLim', [.5 2.5], 'XTick', 1:1:2, 'XTickLabel', {['AS'], ['QS']}); title({[ratIDs2go{ratRUN}] ; ['Resp Peak Interval']}); end

				subplot(nFIGROW, nFIGCOL, 3, 'FontSize', szTXT);
				hold on;
					plot(slRUN, nanstd(sResp(dat2go), [], 1), 'bo');
				hold off;
				if slRUN == 2 set(gca, 'XLim', [.5 2.5], 'XTick', 1:1:2, 'XTickLabel', {['AS'], ['QS']}); title({[ratIDs2go{ratRUN}] ; ['Resp Amplitude variance']}); end

				subplot(nFIGROW, nFIGCOL, 4, 'FontSize', szTXT);
				hold on;
					plot(slRUN, nanmean(zVarMat(dat2go), 2), 'rv');
				hold off;
				if slRUN == 2 set(gca, 'XLim', [.5 2.5], 'XTick', 1:1:2, 'XTickLabel', {['AS'], ['QS']}); title({[ratIDs2go{ratRUN}] ; ['RespB2B variance']}); end

				resp2put = [resp2put [length(dat2go) * SFResp / 1000 / 60 nanmean(diff(find(RespPksIdx(dat2go) == 1))) .* SFResp nanstd(sResp(dat2go), [], 1) nanmean(zVarMat(dat2go), 2)]];
			end %slRUN = 1:1:2

			saveas(picID, [respROOT '\' ratIDs2go{ratRUN} '-RespVar.bmp']); close all;
			respTab = [respTab ; resp2put];

			clear thisRATID preprocLOADROOT preproc slTime sResp tSpace pzlfp m1lfp lfpTSpace midhaAS midhaQS midhaWA midhaASraw midhaQSraw midhaWAraw islerQS islerAS islerQSraw islerASraw resp2SDraw resp2SD pzDelta m1Delta pzDelta1MedRaw pzDelta1Med m1Delta1MedRaw m1Delta1Med pzMUA pzSC pzRasters pzIFR;
		end %ratRUN = 1:1:size(ratIDs2go, 1)
	end %ageRUN = 2:1:size(ageGROUP, 1)
	save([respROOT '\respTab.mat'], 'respTab');
else
	load([respROOT '\respTab.mat']);
end %~exist([respROOT '\respTab.mat'])

respDurationMin = [respTab(respFlag == 1, 1) respTab(respFlag == 1, 5)]; 
respPkInterval = [respTab(respFlag == 1, 2) respTab(respFlag == 1, 6)];
respAmpVar = [respTab(respFlag == 1, 3) respTab(respFlag == 1, 7)];
respIsler = [respTab(respFlag == 1, 4) respTab(respFlag == 1, 8)];

nFIGROW = 1; nFIGCOL = 4;
picID = figure('Color', 'w', 'Position', [25 25 500 450]);
for measRUN = 1:1:nFIGCOL
	if measRUN == 1
		thisDat2plot = respDurationMin; txtTitle2go = ['Resp Peak Interval (n = ' num2str(size(respDurationMin, 1)) ')'];
	elseif measRUN == 2
		thisDat2plot = respPkInterval; txtTitle2go = ['Resp Peak Interval (n = ' num2str(size(respPkInterval, 1)) ')'];
	elseif measRUN == 3
		thisDat2plot = respAmpVar; txtTitle2go = ['RespAmplitude variance (n = ' num2str(size(respAmpVar, 1)) ')'];
	elseif measRUN == 4
		thisDat2plot = respIsler; txtTitle2go = ['RespB2B variance (n = ' num2str(size(respIsler, 1)) ')'];
	end %measRUN == 1

	subplot(nFIGROW, nFIGCOL, measRUN, 'FontSize', szTXT);
	hold on;
		boxplot(thisDat2plot);
		for slRUN = 1:1:2
			for indRUN = 1:1:size(thisDat2plot, 1)
				plot(slRUN, thisDat2plot(indRUN, slRUN), '.', 'Color', ones(1, 3) .* .75);
			end %indRUN = 1:1:size(thisDat2plot, 1)
		end %slRUN = 1:1:2
	hold off;
	set(gca, 'XLim', [.5 2.5], 'XTick', 1:1:2, 'XTickLabel', {['AS'], ['QS']});
	title(txtTitle2go);
end %measRUN = 1:1:2

[h0 p0 ci0 stats0] = ttest2(respDurationMin(:, 1), respDurationMin(:, 2))
[h1 p1 ci1 stats1] = ttest2(respPkInterval(:, 1), respPkInterval(:, 2))
[h2 p2 ci2 stats2] = ttest2(respAmpVar(:, 1), respAmpVar(:, 2))
[h3 p3 ci3 stats3] = ttest2(respIsler(:, 1), respIsler(:, 2))