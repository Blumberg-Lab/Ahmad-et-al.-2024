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
	respROOT = [alzROOT '\RespDeltaPhaPhaV3R1']; if ~exist(respROOT) mkdir(respROOT); end

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

nDegSpace = 30; degSpace = linspace(-pi, pi, nDegSpace);
	indROOT = [respROOT '\nDeg' num2str(nDegSpace)]; if ~exist(indROOT) mkdir(indROOT); end

degSpaceHist = 0:20:360; deg2vis = 10:20:350;

txtOX = {['X'] ; ['O']}; txtLFP = {['PZ'] ; ['M1']};
colLINEs = get(groot, 'DefaultAxesColorOrder'); colLINEs = [colLINEs ; 0 0 1]; colLINEs = [colLINEs ; fliplr(colLINEs)]; szTXT = 7.5;
if ~exist([respROOT '\respPhaPhaHistTab_' num2str(nDegSpace) '.mat'])
	respPhaPhaTab = [];
		respPhaPhaTab1 = [];
		respPhaPhaTab2 = [];
		respPhaPhaTab3 = [];
		respPhaPhaTab4 = [];
	respPhaPhamnTab = [];
		respPhaPhamnTab1 = [];
		respPhaPhamnTab2 = [];
		respPhaPhamnTab3 = [];
		respPhaPhamnTab4 = [];
	respPhaPhaHistTab = [];
		respPhaPhaHistTab1 = [];
		respPhaPhaHistTab2 = [];
		respPhaPhaHistTab3 = [];
		respPhaPhaHistTab4 = [];
	RValTab = [];

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
				pz2HzDelta = jkBandpassFilt(pzlfp, SFtarg, [2 4]);
				pz2HzDeltalsf = downsample(pz2HzDelta, round(length(pz2HzDelta) / length(sResp))); clear pz2HzDelta;
				powPZ2HzLSF = abs(pz2HzDeltalsf) .^ 2; clear pz2HzDeltalsf;
				powPZ2HzLFPIdx = (powPZ2HzLSF > median(powPZ2HzLSF)); clear powPZ2HzLSF;
			m1lfp = preproc.m1lfp;
				m12HzDelta = jkBandpassFilt(m1lfp, SFtarg, [2 4]);
				m12HzDeltalsf = downsample(m12HzDelta, round(length(m12HzDelta) / length(sResp))); clear m12HzDelta;
				powM12HzLSF = abs(m12HzDeltalsf) .^ 2; clear m12HzDeltalsf;
				powM12HzLFPIdx = powM12HzLSF > median(powM12HzLSF); clear powM12HzLSF;
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
				pzlfplsf = downsample(pzDelta, round(length(pzDelta) / length(sResp)));
			m1Delta = preproc.m1Delta;
				m1lfplsf = downsample(m1Delta, round(length(m1Delta) / length(sResp)));
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

			resp2put = []; rvalput = []; rmn2put = []; hist2put = [];
			nFIGROW = 4; nFIGCOL = 6;
			picID = figure('Color', 'w', 'Position', [25 25 1200 720]);
			for slRUN = 1:1:2 	%AS & QS
				if slRUN == 1
					stateDef = [resp2SD ; midhaAS ; islerAS];
					marker2gogo = ['ks'];
				elseif slRUN == 2
					%stateDef = [resp2SD ; midhaQS ; islerQS ; pzDelta1Med ; m1Delta1Med ; powPZ2HzLFPIdx' ; powM12HzLFPIdx'];
					stateDef = [resp2SD ; midhaQS ; islerQS ; pzDelta1Med ; powPZ2HzLFPIdx' ; powM12HzLFPIdx'];
					marker2gogo = ['ro'];
				end %slRUN == 1
				dat2go = find(nansum(stateDef, 1) == size(stateDef, 1));

				phaResp = angle(hilbert(sResp(dat2go)));
				phaPZ = angle(hilbert(pzlfplsf(dat2go)));
				phaM1 = angle(hilbert(m1lfplsf(dat2go)));
					rvalput = [rvalput abs(mean(exp(1i * (phaResp - phaPZ)))) abs(mean(exp(1i * (phaResp - phaM1))))]

				pzrespMap = zeros(nDegSpace, nDegSpace); %zeros(360, 360);
				m1respMap = zeros(nDegSpace, nDegSpace); %zeros(360, 360);

				for datRUN = 1:1:length(phaResp)
					thisRespIdx = dsearchn(transpose(degSpace), phaResp(datRUN));
					thisPZIdx = dsearchn(transpose(degSpace), phaPZ(datRUN));
					thisM1Idx = dsearchn(transpose(degSpace), phaM1(datRUN));

					pzrespMap(thisRespIdx, thisPZIdx) = pzrespMap(thisRespIdx, thisPZIdx) + 1;
					m1respMap(thisRespIdx, thisM1Idx) = m1respMap(thisRespIdx, thisM1Idx) + 1;
					%pzrespMap(thisPZIdx, thisRespIdx) = pzrespMap(thisPZIdx, thisRespIdx) + 1;
					%m1respMap(thisM1Idx, thisRespIdx) = m1respMap(thisM1Idx, thisRespIdx) + 1;

					clear thisRespIdx thisPZIdx thisM1Idx;
				end %datRUN = 1:1:length(phaResp)
				pzrespMap = pzrespMap ./ datRUN;
				m1respMap = m1respMap ./ datRUN;

				subplot(nFIGROW, nFIGCOL, 1 + (slRUN - 1) * nFIGCOL, 'FontSize', szTXT);
				hold on;
					contourf(degSpace, degSpace, pzrespMap, 'LineColor', 'none');
				hold off;
				set(gca, 'XLim', [-pi pi], 'XTick', [-pi pi], 'XTickLabel', {['-pi'], ['pi']}, 'YTick', [-pi pi], 'YTickLabel', {[' '], ['pi']});
				title({[ratIDs2go{ratRUN} ' (' txtOX{1 + pzPLVFlag(ratRUN, 1)} ')'] ; [sleepStage{slRUN} ': Resp-PZ']});

				subplot(nFIGROW, nFIGCOL, 2 + (slRUN - 1) * nFIGCOL, 'FontSize', szTXT);
				hold on;
					contourf(degSpace, degSpace, m1respMap, 'LineColor', 'none');
				hold off;
				set(gca, 'XLim', [-pi pi], 'XTick', [-pi pi], 'XTickLabel', {['-pi'], ['pi']}, 'YTick', [-pi pi], 'YTickLabel', {[' '], ['pi']});
				title({[ratIDs2go{ratRUN} ' (' txtOX{1 + m1PLVFlag(ratRUN, 1)} ')'] ; [sleepStage{slRUN} ': Resp-M1']});

				%Rm:n computation
				%resp & delta ratio == 0.5:1, 1:1, 2:1 (data), 4:1, & 8:1
				xRmn = [.5 1 2 4 8 16];
				rmnMat = nan(2, length(xRmn));
				hist2putMat1 = []; hist2putMat2 = []; rlim1 = []; rlim2 = [];
				for ratioRUN = 1:1:size(rmnMat, 2)
					if ratioRUN == 1 		%.5:1
						dsDivFac = 4;
						thisDSPzDelta = downsample(pzDelta, round(length(pzDelta) / (length(sResp) * dsDivFac)));
							thisDSPzDelta = angle(hilbert(thisDSPzDelta));
							thisDSPzDelta = [thisDSPzDelta(1:length(sResp))];
						thisDSM1Delta = downsample(m1Delta, round(length(m1Delta) / (length(sResp) * dsDivFac)));
							thisDSM1Delta = angle(hilbert(thisDSM1Delta));
							thisDSM1Delta = [thisDSM1Delta(1:length(sResp))];
					elseif ratioRUN == 2 	%1:1
						dsDivFac = 2;
						thisDSPzDelta = downsample(pzDelta, round(length(pzDelta) / (length(sResp) * dsDivFac)));
							thisDSPzDelta = angle(hilbert(thisDSPzDelta));
							thisDSPzDelta = [thisDSPzDelta(1:length(sResp))];
						thisDSM1Delta = downsample(m1Delta, round(length(m1Delta) / (length(sResp) * dsDivFac)));
							thisDSM1Delta = angle(hilbert(thisDSM1Delta));
							thisDSM1Delta = [thisDSM1Delta(1:length(sResp))];
					elseif ratioRUN == 3 	%2:1
						dsDivFac = 1;
						thisDSPzDelta = downsample(pzDelta, round(length(pzDelta) / (length(sResp) * dsDivFac)));
							thisDSPzDelta = angle(hilbert(thisDSPzDelta));
						thisDSM1Delta = downsample(m1Delta, round(length(m1Delta) / (length(sResp) * dsDivFac)));
							thisDSM1Delta = angle(hilbert(thisDSM1Delta));
					elseif ratioRUN == 4 	%4:1
						dsDivFac = .5;
						thisDSPzDelta = downsample(pzDelta, round(length(pzDelta) / (length(sResp) * dsDivFac)));
							thisDSPzDelta = angle(hilbert(thisDSPzDelta));
							thisDSPzDelta = [thisDSPzDelta ; thisDSPzDelta];
						thisDSM1Delta = downsample(m1Delta, round(length(m1Delta) / (length(sResp) * dsDivFac)));
							thisDSM1Delta = angle(hilbert(thisDSM1Delta));
							thisDSM1Delta = [thisDSM1Delta ; thisDSM1Delta];
					elseif ratioRUN == 5
						dsDivFac = .25;
						thisDSPzDelta = downsample(pzDelta, round(length(pzDelta) / (length(sResp) * dsDivFac)));
							thisDSPzDelta = angle(hilbert(thisDSPzDelta));
							thisDSPzDelta = [thisDSPzDelta ; thisDSPzDelta ; thisDSPzDelta ; thisDSPzDelta];
						thisDSM1Delta = downsample(m1Delta, round(length(m1Delta) / (length(sResp) * dsDivFac)));
							thisDSM1Delta = angle(hilbert(thisDSM1Delta));
							thisDSM1Delta = [thisDSM1Delta ; thisDSM1Delta ; thisDSM1Delta ; thisDSM1Delta];
					elseif ratioRUN == 6
						dsDivFac = .125;
						thisDSPzDelta = downsample(pzDelta, round(length(pzDelta) / (length(sResp) * dsDivFac)));
							thisDSPzDelta = angle(hilbert(thisDSPzDelta));
							thisDSPzDelta = [thisDSPzDelta ; thisDSPzDelta ; thisDSPzDelta ; thisDSPzDelta ; thisDSPzDelta ; thisDSPzDelta ; thisDSPzDelta ; thisDSPzDelta];
						thisDSM1Delta = downsample(m1Delta, round(length(m1Delta) / (length(sResp) * dsDivFac)));
							thisDSM1Delta = angle(hilbert(thisDSM1Delta));
							thisDSM1Delta = [thisDSM1Delta ; thisDSM1Delta ; thisDSM1Delta ; thisDSM1Delta ; thisDSM1Delta ; thisDSM1Delta ; thisDSM1Delta ; thisDSM1Delta];
					end %ratioRUN == 1
					thisPZPHAA = thisDSPzDelta(dat2go);
					thisM1PHAA = thisDSM1Delta(dat2go);

					thisHist1 = histc(mod(rad2deg((phaResp - thisPZPHAA)) + 360, 360), degSpaceHist);
						thisHist1(19) = []; thisHist1 = thisHist1 ./ nansum(thisHist1); hist2putMat1 = [hist2putMat1 transpose(thisHist1)];
					thisHist2 = histc(mod(rad2deg((phaResp - thisM1PHAA)) + 360, 360), degSpaceHist);
						thisHist2(19) = []; thisHist2 = thisHist2 ./ nansum(thisHist2); hist2putMat2 = [hist2putMat2 transpose(thisHist2)];

					subplot(nFIGROW, nFIGCOL, ratioRUN + nFIGCOL * 2, 'FontSize', szTXT);
						polarscatter(deg2rad(deg2vis), thisHist1, marker2gogo); hold on;
							rlim1 = [rlim1 ; [nanmin(thisHist1) * .95 nanmax(thisHist1) * 1.05]];
					set(gca, 'rlim', [nanmin(thisHist1) * .95 nanmax(thisHist1) * 1.05]);

					subplot(nFIGROW, nFIGCOL, ratioRUN + nFIGCOL * 3, 'FontSize', szTXT);
						polarscatter(deg2rad(deg2vis), thisHist2, marker2gogo); hold on;
							rlim2 = [rlim2 ; [nanmin(thisHist2) * .95 nanmax(thisHist2) * 1.05]];
					set(gca, 'rlim', [nanmin(thisHist2) * .95 nanmax(thisHist2) * 1.05]);
					
					rmnMat(1, ratioRUN) = abs(mean(exp(1i * (phaResp - thisPZPHAA))));
					rmnMat(2, ratioRUN) = abs(mean(exp(1i * (phaResp - thisM1PHAA))));

					clear dsDivFac thisDSPzDelta thisDSM1Delta thisPZPHAA thisM1PHAA;
				end %ratioRUN = 1:1:size(rmnMat, 2)

				subplot(nFIGROW, nFIGCOL, [3 4] + (slRUN - 1) * nFIGCOL, 'FontSize', szTXT);
				hold on;
					plot(xRmn, rmnMat(1, :));
				hold off;
				set(gca, 'XLim', [.25 max(xRmn) + .5], 'XTick', xRmn); xlabel(['R_m_:_n']); title(['against ' txtLFP{1} ' delta phase']);

				subplot(nFIGROW, nFIGCOL, [5 6] + (slRUN - 1) * nFIGCOL, 'FontSize', szTXT);
				hold on;
					plot(xRmn, rmnMat(2, :));
				hold off;
				set(gca, 'XLim', [.25 max(xRmn) + .5], 'XTick', xRmn); xlabel(['R_m_:_n']); title(['against ' txtLFP{2} ' delta phase']);

				resp2put{slRUN, 1} = pzrespMap; resp2put{slRUN, 2} = m1respMap;
				rmn2put{slRUN, 1} = rmnMat(1, :); rmn2put{slRUN, 2} = rmnMat(2, :);
				hist2put{slRUN, 1} = hist2putMat1; hist2put{slRUN, 2} = hist2putMat2;

				for ratioRUN = 1:1:size(rmnMat, 2)
					subplot(nFIGROW, nFIGCOL, ratioRUN + nFIGCOL * 2, 'FontSize', szTXT);
						set(gca, 'rlim', [nanmin(rlim1(:, 1)) nanmax(rlim1(:, 2))]);
					subplot(nFIGROW, nFIGCOL, ratioRUN + nFIGCOL * 3, 'FontSize', szTXT);
						set(gca, 'rlim', [nanmin(rlim2(:, 1)) nanmax(rlim2(:, 2))]);
				end %ratioRUN = 1:1:size(rmnMat, 2)

				clear stateDef dat2go phaResp phaPZ phaM1 pzrespMap m1respMap;
			end %slRUN = 1:1:2
			savefig([indROOT '\' ratIDs2go{ratRUN} '-RespDeltaPhasePhase.fig']);
			saveas(picID, [indROOT '\' ratIDs2go{ratRUN} '-RespDeltaPhasePhase.bmp']); close all;
			
			respPhaPhaTab1(:, :, ratRUN) = resp2put{1, 1};
			respPhaPhaTab2(:, :, ratRUN) = resp2put{1, 2};
			respPhaPhaTab3(:, :, ratRUN) = resp2put{2, 1};
			respPhaPhaTab4(:, :, ratRUN) = resp2put{2, 2};

			RValTab = [RValTab ; rvalput];

			respPhaPhamnTab1(:, :, ratRUN) = rmn2put{1, 1};
			respPhaPhamnTab2(:, :, ratRUN) = rmn2put{1, 2};
			respPhaPhamnTab3(:, :, ratRUN) = rmn2put{2, 1};
			respPhaPhamnTab4(:, :, ratRUN) = rmn2put{2, 2};

			respPhaPhaHistTab1(:, :, ratRUN) = hist2put{1, 1};
			respPhaPhaHistTab2(:, :, ratRUN) = hist2put{1, 2};
			respPhaPhaHistTab3(:, :, ratRUN) = hist2put{2, 1};
			respPhaPhaHistTab4(:, :, ratRUN) = hist2put{2, 2};

			clear thisRATID preprocLOADROOT preproc slTime sResp tSpace pzlfp m1lfp lfpTSpace midhaAS midhaQS midhaWA midhaASraw midhaQSraw midhaWAraw islerQS islerAS islerQSraw islerASraw resp2SDraw resp2SD pzDelta m1Delta pzDelta1MedRaw pzDelta1Med m1Delta1MedRaw m1Delta1Med pzMUA pzSC pzRasters pzIFR;
		end %ratRUN = 1:1:size(ratIDs2go, 1)
	end %ageRUN = 2:1:size(ageGROUP, 1)
	respPhaPhaTab{1} = respPhaPhaTab1;
	respPhaPhaTab{2} = respPhaPhaTab2; 
	respPhaPhaTab{3} = respPhaPhaTab3; 
	respPhaPhaTab{4} = respPhaPhaTab4;

	respPhaPhamnTab{1} = respPhaPhamnTab1;
	respPhaPhamnTab{2} = respPhaPhamnTab2;
	respPhaPhamnTab{3} = respPhaPhamnTab3;
	respPhaPhamnTab{4} = respPhaPhamnTab4;

	respPhaPhaHistTab{1} = respPhaPhaHistTab1;
	respPhaPhaHistTab{2} = respPhaPhaHistTab2;
	respPhaPhaHistTab{3} = respPhaPhaHistTab3;
	respPhaPhaHistTab{4} = respPhaPhaHistTab4;

	save([respROOT '\respPhaPhaTab_' num2str(nDegSpace) '.mat'], 'respPhaPhaTab');
	save([respROOT '\respPhaPhamnTab_' num2str(nDegSpace) '.mat'], 'respPhaPhamnTab');
	save([respROOT '\respPhaPhaHistTab_' num2str(nDegSpace) '.mat'], 'respPhaPhaHistTab');
	save([respROOT '\RValTab_' num2str(nDegSpace) '.mat'], 'RValTab');
else
	load([respROOT '\respPhaPhaTab_' num2str(nDegSpace) '.mat']);
	load([respROOT '\respPhaPhamnTab_' num2str(nDegSpace) '.mat']);
	load([respROOT '\respPhaPhaHistTab_' num2str(nDegSpace) '.mat']);
	load([respROOT '\RValTab_' num2str(nDegSpace) '.mat']);
end %~exist([respROOT '\respPhaPhaTab.mat'])

xRmn = [.5 1 2 4 8 16];
nFIGROW = 1; nFIGCOL = 2;
figure('Color', 'w');
for slRUN = 2:1:2
	for lfpRUN = 1:1:2
		if lfpRUN == 1
			thisFlag = pzPLVFlag;
		elseif lfpRUN == 2
			thisFlag = m1PLVFlag;
		end %lfpRUN == 1

		thisMap = respPhaPhaTab{lfpRUN + (slRUN - 1) * 2};		

		subplot(nFIGROW, nFIGCOL, lfpRUN, 'FontSize', szTXT);
		hold on;
			contourf(degSpace, degSpace, nanmean(thisMap(:, :, find(thisFlag == 1)), 3), 'LineColor', 'none');
                colormap(gca, parula(64));
		hold off;
		set(gca, 'XLim', [-pi pi], 'XTick', [-pi pi], 'XTickLabel', {['-pi'], ['pi']}, 'YTick', [-pi pi], 'YTickLabel', {[' '], ['pi']});
		if slRUN == 1 & lfpRUN == 1 xlabel(['Delta phase']); ylabel(['respiratory phase']); end
		title([sleepStage{slRUN} ': Resp-' txtLFP{lfpRUN}]);
	end %lfpRUN = 1:1:2
end %slRUN = 2:1:2

nFIGROW = 1; nFIGCOL = 2;
figure('Color', 'w');
for slRUN = 2:1:2
	for lfpRUN = 1:1:2
		if lfpRUN == 1
			thisFlag = pzPLVFlag;
		elseif lfpRUN == 2
			thisFlag = m1PLVFlag;
		end %lfpRUN == 1

		thisRmn = transpose(squeeze(respPhaPhamnTab{lfpRUN + (slRUN - 1) * 2}));
		subplot(nFIGROW, nFIGCOL, lfpRUN, 'FontSize', szTXT);
		hold on;
			plot(xRmn, nanmean(thisRmn(find(thisFlag == 1), :), 1), '-k', 'LineWidth', 1.5);
				plot(xRmn, nanmean(thisRmn(find(thisFlag == 1), :), 1) + nanstd(thisRmn(find(thisFlag == 1), :), [], 1) ./ sqrt(nansum(thisFlag)), ':k', 'LineWidth', .5);
				plot(xRmn, nanmean(thisRmn(find(thisFlag == 1), :), 1) - nanstd(thisRmn(find(thisFlag == 1), :), [], 1) ./ sqrt(nansum(thisFlag)), ':k', 'LineWidth', .5);
		hold off;
		set(gca, 'XLim', [.25 max(xRmn) + .5], 'XTick', xRmn , 'XTickLabel', {['8:1'], ['4:1'], ['2:1'], ['1:1'], ['1:2'], ['1:4']}, 'YLim', [0 .12]); xlabel(['R_m_:_n']); title(['against ' txtLFP{lfpRUN} ' delta phase']);
		clear thisMap thisRmn;
	end %lfpRUN = 1:1:2
end %slRUN = 1:1:2
stophere

nFIGROW = 1; nFIGCOL = 2;
figure('Color', 'w');
for lfpRUN = 1:1:2
	if lfpRUN == 1
		thisFlag = pzPLVFlag;
	elseif lfpRUN == 2
		thisFlag = m1PLVFlag;
	end %lfpRUN == 1

	thisMap = nanmean(respPhaPhaTab{lfpRUN + (2 - 1) * 2}(:, :, find(thisFlag == 1)), 3) - nanmean(respPhaPhaTab{lfpRUN + (1 - 1) * 2}(:, :, find(thisFlag == 1)), 3);
	subplot(nFIGROW, nFIGCOL, lfpRUN, 'FontSize', szTXT);
	hold on;
		contourf(degSpace, degSpace, thisMap, 'LineColor', 'none');
	hold off;
	set(gca, 'XLim', [-pi pi], 'XTick', [-pi pi], 'XTickLabel', {['-pi'], ['pi']}, 'YTick', [-pi pi], 'YTickLabel', {[' '], ['pi']});
	title(['QS - AS: Resp-' txtLFP{lfpRUN}]);
	clear thisMap;
end %lfpRUN = 1:1:2




asRPZ1 = 1; asRPZ50 = 50; asRM11 = 51; asRM150 = 100;
qsRPZ1 = 101; qsRPZ50 = 150; qsRM11 = 151; qsRM150 = 200;

baseFreq = [4 7];
deltaIdx2go = dsearchn(transpose(cohfreq2use), transpose(deltaFreq));
denomIdx2go = dsearchn(transpose(cohfreq2use), transpose(baseFreq));

thisPZMeas = []; thisM1Meas = [];
thisPZMeas2 = []; thisM1Meas2 = []; stdTabPZ = []; stdTabM1 = [];
nFIGROW = 3; nFIGCOL = 2;
picID = figure('Color', 'w', 'Position', [50 50 450 400]);
for slRUN = 1:1:2
	valIdxModi = (slRUN - 1) * 100;
	coh2goRPZ = respCohTab(cohFlag == 1, (asRPZ1 + valIdxModi):(asRPZ50 + valIdxModi));
	coh2goRM1 = respCohTab(cohFlag == 1, (asRM11 + valIdxModi):(asRM150 + valIdxModi));
		for ratRUN = 1:1:size(coh2goRPZ, 1)
			thisMax = max([max(coh2goRPZ(ratRUN, :)) max(coh2goRM1(ratRUN, :))]);
			stdTabPZ(ratRUN, :) = [coh2goRPZ(ratRUN, :) ./ thisMax];
			stdTabM1(ratRUN, :) = [coh2goRM1(ratRUN, :) ./ thisMax];
			clear thisMax;
		end %ratRUN = 1:1:size(coh2goRPZ, 1)

	subplot(nFIGROW, nFIGCOL, [1], 'FontSize', szTXT);
	hold on;
		errorbar(cohfreq2use, nanmean(coh2goRPZ, 1), nanstd(coh2goRPZ, [], 1) ./ sqrt(nansum(cohFlag)));
	hold off;
	set(gca, 'YLim', [0 .025]);
	if slRUN == 2 title([['Resp-PZ (n = ' num2str(nansum(cohFlag)) ')']]); legend(['AS'], ['QS']); end

	subplot(nFIGROW, nFIGCOL, [2], 'FontSize', szTXT);
	hold on;
		errorbar(cohfreq2use, nanmean(coh2goRM1, 1), nanstd(coh2goRM1, [], 1) ./ sqrt(nansum(cohFlag)));
	hold off;
	set(gca, 'YLim', [0 .025]);
	if slRUN == 2 title([['Resp-M1 (n = ' num2str(nansum(cohFlag)) ')']]); legend(['AS'], ['QS']); end

	thisRatioPZ = [];
	thisRatioPZ = nansum(coh2goRPZ(:, deltaIdx2go(1):deltaIdx2go(2)), 2) ./ nansum(coh2goRPZ(:, denomIdx2go(1):denomIdx2go(2)), 2);
	thisPZMeas = [thisPZMeas thisRatioPZ];

	thisRatioM1 = [];
	thisRatioM1 = nansum(coh2goRM1(:, deltaIdx2go(1):deltaIdx2go(2)), 2) ./ nansum(coh2goRM1(:, denomIdx2go(1):denomIdx2go(2)), 2);
	thisM1Meas = [thisM1Meas thisRatioM1];

	thisSTDPZ = [];
	thisSTDPZ = nansum(stdTabPZ(:, deltaIdx2go(1):deltaIdx2go(2)), 2);
	thisPZMeas2 = [thisPZMeas2 thisSTDPZ];

	thisSTDM1 = [];
	thisSTDM1 = nansum(stdTabM1(:, deltaIdx2go(1):deltaIdx2go(2)), 2);
	thisM1Meas2 = [thisM1Meas2 thisSTDM1];

	for typeRUN = 1:1:2
		for measRUN = 1:1:2
			if typeRUN == 1
				ymeasLab = ['delta / theta ratio'];
				if measRUN == 1
					ttMeas = thisPZMeas;
				elseif measRUN == 2
					ttMeas = thisM1Meas;
				end %measRUN == 1
			else
				ymeasLab = ['std. delta coh'];
				if measRUN == 1
					ttMeas = thisPZMeas2;
				elseif measRUN == 2
					ttMeas = thisM1Meas2;
				end %measRUN == 1
			end %typeRUN == 1
			subplot(nFIGROW, nFIGCOL, [3 + (measRUN - 1) + (typeRUN - 1) * nFIGCOL], 'FontSize', szTXT);
			hold on;
				if slRUN == 2
					boxplot(ttMeas);

					for slRUN = 1:1:2
						for indRUN = 1:1:size(ttMeas, 1)
							plot(slRUN, ttMeas(indRUN, slRUN), '.', 'Color', ones(1, 3) .* .75);
						end %indRUN = 1:1:size(ttMeas, 1)
					end %slRUN = 1:1:2
					set(gca, 'XLim', [.5 2.5], 'XTick', 1:1:2, 'XTickLabel', {['AS'], ['QS']}); ylabel(ymeasLab);
				end %slRUN == 2
			hold off;
		end %measRUN = 1:1:2
	end %typeRUN = 1:1:2
end %slRUN = 1:1:2

[h1 p1 ci1 stats1] = ttest2(thisPZMeas(:, 1), thisPZMeas(:, 2))
[h2 p2 ci2 stats2] = ttest2(thisM1Meas(:, 1), thisM1Meas(:, 2))
[h3 p3 ci3 stats3] = ttest2(thisPZMeas2(:, 1), thisPZMeas2(:, 2))
[h4 p4 ci4 stats4] = ttest2(thisM1Meas2(:, 1), thisM1Meas2(:, 2))
here

respPkInterval = [respTab(cohFlag == 1, 1) respTab(cohFlag == 1, 3)];
respAmpVar = [respTab(cohFlag == 1, 2) respTab(cohFlag == 1, 4)];

nFIGROW = 1; nFIGCOL = 2;
picID = figure('Color', 'w', 'Position', [25 25 500 450]);
for measRUN = 1:1:2
	if measRUN == 1
		thisDat2plot = respPkInterval; txtTitle2go = ['Resp Peak Interval (n = ' num2str(size(respPkInterval, 1)) ')'];
	elseif measRUN == 2
		thisDat2plot = respAmpVar; txtTitle2go = ['RespAmplitude variance (n = ' num2str(size(respAmpVar, 1)) ')'];
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

[h1 p1 ci1 stats1] = ttest2(respPkInterval(:, 1), respPkInterval(:, 2))
[h2 p2 ci2 stats2] = ttest2(respAmpVar(:, 1), respAmpVar(:, 2))