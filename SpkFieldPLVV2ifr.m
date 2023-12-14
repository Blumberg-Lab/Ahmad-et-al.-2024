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
	plvROOT = [alzROOT '\spkfieldV2ifr']; if ~exist(plvROOT) mkdir(plvROOT); end
		visROOT = [plvROOT '\VisualSummary']; if ~exist(visROOT) mkdir(visROOT); end

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

nITER = 1000; degSpace = 0:20:360;

txtOX = {['X'] ; ['O']}; txtLFPS = {['PZ'] ; ['M1']}; txtSIG = {[' '] ; ['*']};
colLINEs = get(groot, 'DefaultAxesColorOrder'); colLINEs = [colLINEs ; 0 0 1]; colLINEs = [colLINEs ; fliplr(colLINEs)]; szTXT = 7.5;
if ~exist([plvROOT '\plvMasterTab.mat'])
	plvMasterTab = []; %50 x 2; first 50 for AS & secon 50 for QS

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
				smpzIFR = smooth(pzIFR);
				envpzIFR = envelope(smpzIFR, 500, 'analytic');
				smpzIFRAbvThr = nan(length(smpzIFR), 1);
				smpzIFRAbvThr(find(smpzIFR > (median(envpzIFR) * 1))) = smpzIFR(find(smpzIFR > (median(envpzIFR) * 1)));

			nFIGROW = 4; nFIGCOL = 3;
			plv2put = []; %6 stats + 18 histogram values * 2 states * 2 lfps * 3 cris
			picID = figure('Color', 'w', 'Position', [50 50 1470 800], 'Visible', 1);
			for slRUN = 1:1:2 %AS & QS
				if slRUN == 1
					stateDef = [resp2SDraw ; midhaASraw ; islerASraw];
				elseif slRUN == 2
					%stateDef = [resp2SDraw ; midhaQSraw ; islerQSraw ; pzDelta1MedRaw ; m1Delta1MedRaw];
                    stateDef = [resp2SDraw ; midhaQSraw ; islerQSraw ; pzDelta1MedRaw];
				end %slRUN == 1

				for lfpRUN = 1:1:2
					if lfpRUN == 1
						thisDelta = pzDelta;
						thisFlag = pzPLVFlag;
					elseif lfpRUN == 2
						thisDelta = m1Delta;
						thisFlag = m1PLVFlag;
					end %lfpRUN == 1

					thisDeltaPhaseRadian = angle(hilbert(thisDelta));	%in radian
					thisDeltaPhaseDeg = mod(rad2deg(thisDeltaPhaseRadian) + 360, 360);

					for criRUN = 1:1:nFIGCOL
						pzIFRidx = zeros(1, length(resp2SDraw));
						if criRUN == 1
							pzIFRidx(find(smpzIFRAbvThr > 0)) = 1;
						elseif criRUN == 2
							pzIFRidx(find(smpzIFRAbvThr > .25)) = 1;
						elseif criRUN == 3
							pzIFRidx(find(smpzIFRAbvThr > .50)) = 1;
						elseif criRUN == 4
							pzIFRidx(find(smpzIFRAbvThr > .75)) = 1;
						elseif criRUN == 5
							pzIFRidx(find(smpzIFRAbvThr > .90)) = 1;
						end %criRUN == 1

						stateDefhere = [stateDef ; pzIFRidx];
						dat2go = find(nansum(stateDefhere, 1) == size(stateDefhere, 1));

						[thisZPLV thisRawPLV] = jkGetZPLV(dat2go, thisDeltaPhaseDeg, nITER);

						if length(dat2go) == 0
							rayleighP = nan; rayleighZ = nan;
							meanPhaseRadian = nan; meanKappa = nan;
							thisHist2go = nan(1, length(degSpace));
						else
							[rayleighP rayleighZ] = circ_rtest(thisDeltaPhaseRadian(dat2go)); if rayleighP < .05 thisSIG = txtSIG{2}; else thisSIG = txtSIG{1}; end
							meanPhaseRadian = circ_mean(thisDeltaPhaseRadian(dat2go));
							meanKappa = circ_kappa(thisDeltaPhaseRadian(dat2go));

							thisHist2go = histc(thisDeltaPhaseDeg(dat2go), degSpace) ./ nansum(histc(thisDeltaPhaseDeg(dat2go), degSpace));
							if size(thisHist2go, 1) > size(thisHist2go, 2)
								thisHist2go = transpose(thisHist2go);
							end %size(thisHist2go, 1) > size(thisHist2go, 2)
						end %length(dat2go) == 0
						thisHist2go(end) = [];

						plv2put = [plv2put ... 
									[thisRawPLV thisZPLV rayleighP rayleighZ meanPhaseRadian meanKappa thisHist2go]];

						subplot(nFIGROW, nFIGCOL, criRUN + (slRUN - 1) * (nFIGCOL * 2) + (lfpRUN - 1) * nFIGCOL, 'FontSize', szTXT);
						hold on;
							bar(10:20:350, thisHist2go);
						hold off;
						set(gca, 'XLim', [-.5 360.5], 'XTick', [0 180 360]);
						xlabel(['Phase (deg)']); if lfpRUN == 1 ylabel(['observations (%)']); end
						title({[thisRATID ' (' txtOX{1 + thisFlag(ratRUN)} ')'] ; [sleepStage{slRUN} '-' txtLFPS{lfpRUN} ' delta'] ; ... 
								[thisSIG ': raw PLV = ' df2str(thisRawPLV) '; z-PLV = ' df2str(thisZPLV) '; mean Rad = ' df2str(mod(rad2deg(meanPhaseRadian) + 360, 360))]});

						clear stateDefhere dat2go 
					end %criRUN = 1:1:5
					clear thisDelta thisDeltaPhaseRadian thisDeltaPhaseDeg thisZPLV thisRawPLV rayleighP rayleighZ meanPhaseRadian meanKappa thisHist2go;
				end %lfpRUN = 1:1:2
			end %slRUN = 1:1:2 %AS & QS
			saveas(picID, [visROOT '\' thisRATID '-PZIfr-PLV.bmp']); close all;
			plvMasterTab = [plvMasterTab ; [ratRUN plv2put]];

			clear thisRATID preprocLOADROOT preproc slTime sResp tSpace pzlfp m1lfp lfpTSpace midhaAS midhaQS midhaWA midhaASraw midhaQSraw midhaWAraw islerQS islerAS islerQSraw islerASraw resp2SDraw resp2SD pzDelta m1Delta pzDelta1MedRaw pzDelta1Med m1Delta1MedRaw m1Delta1Med pzMUA pzSC pzRasters pzIFR;
		end %ratRUN = 1:1:size(ratIDs2go, 1)
	end %ageRUN = 2:1:size(ageGROUP, 1)
	save([plvROOT '\plvMasterTab.mat'], 'plvMasterTab');
else
	load([plvROOT '\plvMasterTab.mat']);
end %~exist([plvROOT '\plvMasterTab.mat'])

pRATID = 1; pRAWPLV_ASPZ = 2; pZPLV_ASPZ = 3; pPVAL_ASPZ = 4; pMeanPhase_ASPZ = 6; pMeanKappa_ASPZ = 7; pHIST1_ASPZ = 8; pHIST18_ASPZ = 25;
				%26 ~ 49 for ASPZ cri > .25
				%50 ~ 73 for ASPZ cri > .5
			pRAWPLV_ASM1 = 74; pZPLV_ASM1 = 75; pPVAL_ASM1 = 76; pMeanPhase_ASM1 = 77; pMeanKappa_ASM1 = 78; pHIST1_ASM1 = 79; pHIST18_ASM1 = 97;
				%98 ~ 121 for ASM1 cri > .25
				%122 ~ 145 for ASM1 cri > .5
			pRAWPLV_QSPZ = 146; pZPLV_QSPZ = 147; pPVAL_QSPZ = 148; pMeanPhase_QSPZ = 149; pMeanKappa_QSPZ = 150; pHIST1_QSPZ = 151; pHIST18_QSPZ = 169;
				%170 ~ 193 for QSPZ cri > .25
				%194 ~ 217 for QSPZ cri > .5
			pRAWPLV_QSM1 = 218; pZPLV_QSM1 = 219; pPVAL_QSM1 = 220; pMeanPhase_QSM1 = 221; pMeanKappa_QSM1 = 222; pHIST1_QSM1 = 223; pHIST18_QSM1 = 241;
				%242 ~ 265 for QSM1 cri > .25
				%266 ~ 289 for QSM1 cri > .5

nFIGROW = 2; nFIGCOL = 4;
%population [visual check]
for criRUN = 3%1:1:3
	deg2plot = 10:20:350; boxplotTab = []; boxplotTab2 = []; cnt = 1;
	picID = figure('Color', 'w', 'Position', [50 50 1150 530]);
	for slRUN = 1:1:2
		for lfpRUN = 1:1:2
			if lfpRUN == 1
				thisPLVFlag = find(pzPLVFlag == 1);
	 		elseif lfpRUN == 2
	 			thisPLVFlag = find(m1PLVFlag == 1);
			end %lfpRUN == 1
			thisModi = (lfpRUN - 1) * 72 + (slRUN - 1) * 144 + (criRUN - 1) * 24; nRats2use = length(thisPLVFlag);

			subplot(nFIGROW, nFIGCOL, [1 2], 'FontSize', szTXT);
			hold on;
				thisChk = plvMasterTab(thisPLVFlag, pZPLV_ASPZ + thisModi); boxplotTab = [boxplotTab ; [thisChk ones(length(thisChk), 1) .* cnt]];
				for indRUN = 1:1:size(thisChk, 1)
					plot(lfpRUN + (slRUN - 1) * 2, thisChk(indRUN, 1), '.', 'Color', ones(1, 3) .* .75);
				end %indRUN = 1:1:size(thisChk, 1)
				plot(lfpRUN + (slRUN - 1) * 2, nanmean(thisChk, 1), 'ks');
				[hhh ppp] = ttest(thisChk, 1.96)
				%if lfpRUN == 2 & slRUN == 2 boxplot(boxplotTab(:, 1), boxplotTab(:, 2)); end
			hold off;
			set(gca, 'XLim', [.5 4.5], 'XTick', 1:1:4, 'XTickLabel', {['PZ-AS'], ['M1-AS'], ['PZ-QS'], ['M1-QS']});%, 'YLim', [-2.5 15]);
			% ylabel({['Mean raw PLV']});
			ylabel({['Mean Z-PLV']});
			title({['Population [PZ: ' num2str(nansum(pzPLVFlag, 1)) '; M1: ' num2str(nansum(m1PLVFlag, 1)) ']']});

			subplot(nFIGROW, nFIGCOL, 3, 'FontSize', szTXT);
			hold on;
				thisChk2 = plvMasterTab(thisPLVFlag, pMeanPhase_ASPZ + thisModi);
				sum(isnan(thisChk2))
				thisChk2(isnan(thisChk2)) = [];
				errorbar(lfpRUN + (slRUN - 1) * 2, mod(rad2deg(circ_mean((thisChk2))) + 720, 360), ... 
				 									mod(rad2deg(circ_std((thisChk2))) + 720, 360) ./ sqrt(nRats2use), 'ks');
			hold off;
			set(gca, 'XLim', [.5 4.5], 'XTick', 1:1:4, 'XTickLabel', {['PZ_A_S'], ['M1_A_S'], ['PZ_Q_S'], ['M1_Q_S']}, 'YLim', [-5 365], 'YTick', 0:180:360, 'YTickLabel', {['Peak'], ['Trough'], ['Peak']});
			ylabel({['mean phase (deg)']}); if lfpRUN == 1 & slRUN == 1 title(['NOTE: data are circular!']); end

			subplot(nFIGROW, nFIGCOL, 4, 'FontSize', szTXT);
			hold on;
				thisChk3 = plvMasterTab(thisPLVFlag, pMeanKappa_ASPZ + thisModi); boxplotTab2 = [boxplotTab2 ; [thisChk3 ones(length(thisChk3), 1) .* cnt]]; cnt = cnt + 1;
				for indRUN = 1:1:size(thisChk3, 1)
					plot(lfpRUN + (slRUN - 1) * 2, thisChk3(indRUN, 1), '.', 'Color', ones(1, 3) .* .75);
				end %indRUN = 1:1:size(thisChk2, 1)
				plot(lfpRUN + (slRUN - 1) * 2, nanmean(thisChk3, 1), 'ks');
				%if lfpRUN == 2 & slRUN == 2 boxplot(boxplotTab2(:, 1), boxplotTab2(:, 2)); end
			hold off;
			set(gca, 'XLim', [.5 4.5], 'XTick', 1:1:4, 'XTickLabel', {['PZ-AS'], ['M1-AS'], ['PZ-QS'], ['M1-QS']}); %, 'YLim', [0 1]);
			ylabel({['mean kappa'] ; ['(phase-locking concentration)']});

			subplot(nFIGROW, nFIGCOL, nFIGCOL + lfpRUN + (slRUN - 1) * 2, 'FontSize', szTXT);
			hold on;
				plot(deg2plot, nanmean(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), 1), 'k-', 'LineWidth', 2);
					plot(deg2plot, nanmean(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), 1) + nanstd(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), [], 1) / sqrt(nRats2use), 'k:', 'LineWidth', 1);
					plot(deg2plot, nanmean(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), 1) - nanstd(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), [], 1) / sqrt(nRats2use), 'k:', 'LineWidth', 1);
			hold off;
			set(gca, 'XLim', [5 355], 'XTick', [180], 'YLim', [0 .2]); %, 'YLim', [0 2.5]);
			xlabel(['Phase (deg)']); if lfpRUN == 1 & slRUN == 1 ylabel(['AVG observations (%)']); end
			title([sleepStage{slRUN} '-' txtLFPS{lfpRUN}]);
		end %lfpRUN = 1:1:2
	end %slRUN = 1:1:2
end %criRUN = 1:1:3

% nFIGROW = 2; nFIGCOL = 4;
% %population [visual check]
% for criRUN = 3%1:1:3
% 	deg2plot = 10:20:350; boxplotTab = []; boxplotTab2 = []; cnt = 1;
% 	picID = figure('Color', 'w', 'Position', [50 50 1150 530]);
% 	for slRUN = 1:1:2
% 		for lfpRUN = 1:1:2
% 			if lfpRUN == 1
% 				thisPLVFlag = find(pzPLVFlag == 1);
% 	 		elseif lfpRUN == 2
% 	 			thisPLVFlag = find(m1PLVFlag == 1);
% 			end %lfpRUN == 1
% 			thisModi = (lfpRUN - 1) * 72 + (slRUN - 1) * 144 + (criRUN - 1) * 24; nRats2use = length(thisPLVFlag);

% 			subplot(nFIGROW, nFIGCOL, [1 2], 'FontSize', szTXT);
% 			hold on;
% 				thisChk = plvMasterTab(thisPLVFlag, pZPLV_ASPZ + thisModi); boxplotTab = [boxplotTab ; [thisChk ones(length(thisChk), 1) .* cnt]];
% 				for indRUN = 1:1:size(thisChk, 1)
% 					plot(lfpRUN + (slRUN - 1) * 2, thisChk(indRUN, 1), '.', 'Color', ones(1, 3) .* .75);
% 				end %indRUN = 1:1:size(thisChk, 1)
% 				plot(lfpRUN + (slRUN - 1) * 2, nanmean(thisChk, 1), 'ks');
% 				if lfpRUN == 2 & slRUN == 2 boxplot(boxplotTab(:, 1), boxplotTab(:, 2)); end
% 			hold off;
% 			set(gca, 'XLim', [.5 4.5], 'XTick', 1:1:4, 'XTickLabel', {['PZ-AS'], ['M1-AS'], ['PZ-QS'], ['M1-QS']});%, 'YLim', [-2.5 15]);
% 			% ylabel({['Mean raw PLV']});
% 			ylabel({['Mean Z-PLV']});
% 			title({['Population [PZ: ' num2str(nansum(pzPLVFlag, 1)) '; M1: ' num2str(nansum(m1PLVFlag, 1)) ']']});

% 			subplot(nFIGROW, nFIGCOL, 3, 'FontSize', szTXT);
% 			hold on;
% 				thisChk2 = plvMasterTab(thisPLVFlag, pMeanPhase_ASPZ + thisModi);
% 				thisChk2(isnan(thisChk2)) = [];
% 				errorbar(lfpRUN + (slRUN - 1) * 2, mod(rad2deg(circ_mean((thisChk2))) + 720, 360), ... 
% 				 									mod(rad2deg(circ_std((thisChk2))) + 720, 360) ./ sqrt(nRats2use), 'ks');
% 			hold off;
% 			set(gca, 'XLim', [.5 4.5], 'XTick', 1:1:4, 'XTickLabel', {['PZ_A_S'], ['M1_A_S'], ['PZ_Q_S'], ['M1_Q_S']}, 'YLim', [-5 365], 'YTick', 0:180:360, 'YTickLabel', {['Peak'], ['Trough'], ['Peak']});
% 			ylabel({['mean phase (deg)']}); if lfpRUN == 1 & slRUN == 1 title(['NOTE: data are circular!']); end

% 			subplot(nFIGROW, nFIGCOL, 4, 'FontSize', szTXT);
% 			hold on;
% 				thisChk3 = plvMasterTab(thisPLVFlag, pMeanKappa_ASPZ + thisModi); boxplotTab2 = [boxplotTab2 ; [thisChk3 ones(length(thisChk3), 1) .* cnt]]; cnt = cnt + 1;
% 				for indRUN = 1:1:size(thisChk3, 1)
% 					plot(lfpRUN + (slRUN - 1) * 2, thisChk3(indRUN, 1), '.', 'Color', ones(1, 3) .* .75);
% 				end %indRUN = 1:1:size(thisChk2, 1)
% 				plot(lfpRUN + (slRUN - 1) * 2, nanmean(thisChk3, 1), 'ks');
% 				if lfpRUN == 2 & slRUN == 2 boxplot(boxplotTab2(:, 1), boxplotTab2(:, 2)); end
% 			hold off;
% 			set(gca, 'XLim', [.5 4.5], 'XTick', 1:1:4, 'XTickLabel', {['PZ-AS'], ['M1-AS'], ['PZ-QS'], ['M1-QS']}); %, 'YLim', [0 1]);
% 			ylabel({['mean kappa'] ; ['(phase-locking concentration)']});

% 			subplot(nFIGROW, nFIGCOL, nFIGCOL + lfpRUN + (slRUN - 1) * 2, 'FontSize', szTXT);
% 			hold on;
% 				plot(deg2plot, nanmean(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), 1), 'k-', 'LineWidth', 2);
% 					plot(deg2plot, nanmean(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), 1) + nanstd(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), [], 1) / sqrt(nRats2use), 'k:', 'LineWidth', 1);
% 					plot(deg2plot, nanmean(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), 1) - nanstd(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), [], 1) / sqrt(nRats2use), 'k:', 'LineWidth', 1);
% 			hold off;
% 			set(gca, 'XLim', [5 355], 'XTick', [180]); %, 'YLim', [0 2.5]);
% 			xlabel(['Phase (deg)']); if lfpRUN == 1 & slRUN == 1 ylabel(['AVG observations (%)']); end
% 			title([sleepStage{slRUN} '-' txtLFPS{lfpRUN}]);
% 		end %lfpRUN = 1:1:2
% 	end %slRUN = 1:1:2
% end %criRUN = 1:1:3
fafsafas
%decided to visualize QS, > .5 cri results only
criRUN = 3; slRUN = 2;
nFIGROW = 2; nFIGCOL = 2;
deg2plot = 10:20:350; boxplotTab = []; boxplotTab2 = []; cnt = 1;
picID = figure('Color', 'w', 'Position', [50 50 750 530]);
for lfpRUN = 1:1:2
	if lfpRUN == 1
		thisPLVFlag = find(pzPLVFlag == 1);
	elseif lfpRUN == 2
		thisPLVFlag = find(m1PLVFlag == 1);
	end %lfpRUN == 1
	thisModi = (lfpRUN - 1) * 72 + (slRUN - 1) * 144 + (criRUN - 1) * 24; nRats2use = length(thisPLVFlag);

	subplot(nFIGROW, nFIGCOL, [1], 'FontSize', szTXT);
	hold on;
		thisChk = plvMasterTab(thisPLVFlag, pZPLV_ASPZ + thisModi); boxplotTab = [boxplotTab ; [thisChk ones(length(thisChk), 1) .* cnt]];
		for indRUN = 1:1:size(thisChk, 1)
			plot(lfpRUN, thisChk(indRUN, 1), '.', 'Color', ones(1, 3) .* .75);
		end %indRUN = 1:1:size(thisChk, 1)
		if lfpRUN == 2 & slRUN == 2 boxplot(boxplotTab(:, 1), boxplotTab(:, 2)); end
	hold off;
	set(gca, 'XLim', [.5 2.5], 'XTick', 1:1:2, 'XTickLabel', {['PZ-QS'], ['M1-QS']});%, 'YLim', [-2.5 15]);
	% ylabel({['Mean raw PLV']});
	ylabel({['Mean Z-PLV']});
	title({['Population [PZ: ' num2str(nansum(pzPLVFlag, 1)) '; M1: ' num2str(nansum(m1PLVFlag, 1)) ']']});

	subplot(nFIGROW, nFIGCOL, [2], 'FontSize', szTXT);
	hold on;
		thisChk2 = plvMasterTab(thisPLVFlag, pMeanKappa_ASPZ + thisModi); boxplotTab2 = [boxplotTab2 ; [thisChk2 ones(length(thisChk2), 1) .* cnt]]; cnt = cnt + 1;
		for indRUN = 1:1:size(thisChk2, 1)
			plot(lfpRUN, thisChk2(indRUN, 1), '.', 'Color', ones(1, 3) .* .75);
		end %indRUN = 1:1:size(thisChk2, 1)
		if lfpRUN == 2 & slRUN == 2 boxplot(boxplotTab2(:, 1), boxplotTab2(:, 2)); end
	hold off;
	set(gca, 'XLim', [.5 2.5], 'XTick', 3:1:4, 'XTickLabel', {['PZ-QS'], ['M1-QS']}); %, 'YLim', [0 1]);
	ylabel({['mean kappa'] ; ['(phase-locking concentration)']});

	subplot(nFIGROW, nFIGCOL, lfpRUN + nFIGCOL, 'FontSize', szTXT);
	hold on;
		plot(deg2plot, nanmean(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), 1), 'k-', 'LineWidth', 2);
			plot(deg2plot, nanmean(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), 1) + nanstd(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), [], 1) / sqrt(nRats2use), 'k:', 'LineWidth', 1);
			plot(deg2plot, nanmean(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), 1) - nanstd(plvMasterTab(thisPLVFlag, ((pHIST1_ASPZ:pHIST18_ASPZ) + thisModi)), [], 1) / sqrt(nRats2use), 'k:', 'LineWidth', 1);
	hold off;
	set(gca, 'XLim', [5 355], 'XTick', [180]); %, 'YLim', [0 2.5]);
	xlabel(['Phase (deg)']); if lfpRUN == 1 & slRUN == 1 ylabel(['AVG observations (%)']); end
	title([sleepStage{slRUN} '-' txtLFPS{lfpRUN}]);
end %lfpRUN = 1:1:2
sssj
saveas(picID, [plvROOT '\' 'Population-PLVSummaryV2.bmp']);
print -depsc -tiff -r300 -painters 'Pop-PZIFR-Delta-PLV.eps';

%population
picID = figure('Color', 'w', 'Position', [50 50 1150 530]);
for slRUN = 1:1:2
	for lfpRUN = 1:1:2
		if lfpRUN == 1
			thisPLVFlag = find(pzPLVFlag == 1);
 		elseif lfpRUN == 2
 			thisPLVFlag = find(m1PLVFlag == 1);
		end %lfpRUN == 1
		thisModi = (lfpRUN - 1) * 25 + (slRUN - 1) * 50; nRats2use = length(thisPLVFlag);

		subplot(nFIGROW, nFIGCOL, 1, 'FontSize', szTXT);
		hold on;
			yyaxis left;
				errorbar(lfpRUN + (slRUN - 1) * 2, nanmean(ratTab(thisPLVFlag, rRPLV_ASPZ + thisModi), 1), ... 
													nanstd(ratTab(thisPLVFlag, rRPLV_ASPZ + thisModi), [], 1) ./ sqrt(nRats2use), 'bs');
			yyaxis right;
				errorbar(lfpRUN + (slRUN - 1) * 2, nanmean(ratTab(thisPLVFlag, rZPLV_ASPZ + thisModi), 1), ... 
													nanstd(ratTab(thisPLVFlag, rZPLV_ASPZ + thisModi), [], 1) ./ sqrt(nRats2use), 'ro');
				if slRUN == 1 & lfpRUN == 1 ylabel({['Mean Z-PLV']}); end
			yyaxis left;
		hold off;
		set(gca, 'XLim', [.5 4.5], 'XTick', 1:1:4, 'XTickLabel', {['PZ_A_S'], ['M1_A_S'], ['PZ_Q_S'], ['M1_Q_S']});
		ylabel({['Mean raw PLV']});
		ylabel({['Mean Z-PLV']});
		title({['Population [PZ: ' num2str(nansum(pzPLVFlag, 1)) '; M1: ' num2str(nansum(m1PLVFlag, 1)) ']']});

		subplot(nFIGROW, nFIGCOL, 2, 'FontSize', szTXT);
		hold on;
			yyaxis left;
				errorbar(lfpRUN + (slRUN - 1) * 2, nanmean(ratTab(thisPLVFlag, (rSIGCELLS_ASPZ + thisModi)), 1), ... 
													nanstd(ratTab(thisPLVFlag, (rSIGCELLS_ASPZ + thisModi)), 1) ./ sqrt(nRats2use), 'bs');
			yyaxis right;
				errorbar(lfpRUN + (slRUN - 1) * 2, (nanmean(ratTab(thisPLVFlag, (rSIGPERC_ASPZ + thisModi)), 1)), ... 
													nanstd(ratTab(thisPLVFlag, (rSIGPERC_ASPZ + thisModi)), 1) ./ sqrt(nRats2use), 'ro');
				if slRUN == 1 & lfpRUN == 1 ylabel({['Prop sig cells']}); end
			yyaxis left;
		hold off;
		set(gca, 'XLim', [.5 4.5], 'XTick', 1:1:4, 'XTickLabel', {['PZ_A_S'], ['M1_A_S'], ['PZ_Q_S'], ['M1_Q_S']});
		ylabel({['Avg. # of sig cells']});

		subplot(nFIGROW, nFIGCOL, 3, 'FontSize', szTXT);
		hold on;
			errorbar(lfpRUN + (slRUN - 1) * 2, mod(rad2deg(circ_mean(deg2rad(ratTab(thisPLVFlag, rMeanPhase_ASPZ + thisModi)))) + 720, 360), ... 
			 									mod(rad2deg(circ_std(deg2rad(ratTab(thisPLVFlag, rMeanPhase_ASPZ + thisModi)))) + 720, 360) ./ sqrt(nRats2use), 'ks');
		hold off;
		set(gca, 'XLim', [.5 4.5], 'XTick', 1:1:4, 'XTickLabel', {['PZ_A_S'], ['M1_A_S'], ['PZ_Q_S'], ['M1_Q_S']}, 'YLim', [-5 365], 'YTick', 0:180:360, 'YTickLabel', {['Peak'], ['Trough'], ['Peak']});
		ylabel({['mean phase (deg)']}); if lfpRUN == 1 & slRUN == 1 title(['NOTE: data are circular!']); end

		subplot(nFIGROW, nFIGCOL, 4, 'FontSize', szTXT);
		hold on;
			errorbar(lfpRUN + (slRUN - 1) * 2, nanmean(ratTab(thisPLVFlag, rMeanKappa_ASPZ + thisModi), 1), ... 
												nanstd(ratTab(thisPLVFlag, rMeanKappa_ASPZ + thisModi), [], 1) ./ sqrt(nRats2use), 'ks');
		hold off;
		set(gca, 'XLim', [.5 4.5], 'XTick', 1:1:4, 'XTickLabel', {['PZ_A_S'], ['M1_A_S'], ['PZ_Q_S'], ['M1_Q_S']});
		ylabel({['mean kappa'] ; ['(phase-locking concentration)']});

		subplot(nFIGROW, nFIGCOL, nFIGCOL + lfpRUN + (slRUN - 1) * 2, 'FontSize', szTXT);
		hold on;
			errorbar(degSpace, nanmean(ratTab(thisPLVFlag, ((rHIST1_ASPZ:rHIST19_ASPZ) + thisModi)), 1), ... 
								nanstd(ratTab(thisPLVFlag, ((rHIST1_ASPZ:rHIST19_ASPZ) + thisModi)), [], 1) / sqrt(nRats2use), '-b');
		hold off;
		set(gca, 'XLim', [-15 355], 'XTick', [0:180:360]);
		xlabel(['Phase (deg)']); if lfpRUN == 1 & slRUN == 1 ylabel(['AVG observations (%)']); end
		title([sleepStage{slRUN} '-' txtLFPS{lfpRUN}]);
	end %lfpRUN = 1:1:2
end %slRUN = 1:1:2
saveas(picID, [plvROOT '\' 'Population-PLVSummary.bmp']);
print -depsc -tiff -r300 -painters 'Pop-PLV.eps';