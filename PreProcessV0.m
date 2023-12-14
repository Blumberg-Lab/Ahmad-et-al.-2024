%Sleep project by Blumberg lab
%Midha, Greta, & Mark
%
%Jangjin Kim, 2023-Aug-16/2023-Oct-26 [for P10 data]

%initialization
clear all; close all; fclose all; clc;

%def basic paths
datROOT = ['\\lc-rs-store21.hpc.uiowa.edu\Blumberg_Lab_LSS\Jin\Data\'];
alzROOT = ['G:\Blumberg\PZProjectV3']; if ~exist(alzROOT) mkdir(alzROOT); end
	preprocROOT = [alzROOT '\preprocessV0']; if ~exist(preprocROOT) mkdir(preprocROOT); end
		visSum = [preprocROOT '\VisualSumamry']; if ~exist(visSum) mkdir(visSum); end

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
thetaFreq = [4 10];

if ~exist([visSum '\qsPropTab.mat'])
	qsPropTab = [];
	colLINEs = get(groot, 'DefaultAxesColorOrder'); colLINEs = [colLINEs ; 0 0 1]; colLINEs = [colLINEs ; fliplr(colLINEs)]; szTXT = 7.5;
	for ageRUN = 1:1:1%size(ageGROUP, 1)
		ratDat = dir([datROOT '\' ageGROUP{ageRUN} '\' ageGROUP{ageRUN} '*.mat']);
		for ratRUN = 1%1:1:size(ratDat, 1)	%size(ratDat, 1):-1:1
			whereUnderscore = findstr(ratDat(ratRUN).name, '_');
			thisRATID = ratDat(ratRUN).name(whereUnderscore+1:end-4); clear whereUnderscore
			disp(['Working on ' thisRATID])

			%channels to use for analyses [info. from Midha]
			if ageRUN == 1
				switch thisRATID
					case ['10D1']		%10D1
						pz2go = [16]; %from 6-16
						m12go = [3]; %from 2-16

						SFclfp = 1017.3;
					case ['12D2'] 		%12D2
						pz2go = [9]; %from 1-9
						m12go = [8]; %from 8-16

						SFclfp = 1017.3;
					case ['3J1']
						pz2go = [16]; %from 2-16
						m12go = [5]; %from 4-16
					case ['6J3']
						pz2go = [15]; %from 9, 10, 11, & 15
						m12go = [8]; %from 8, 10, 12-16
					case ['7B1']
						pz2go = [15]; %from 4, 6-15
						m12go = [11]; %from 11-16
					case ['7F1']
						pz2go = [10]; %from 2-4, 6-10
						m12go = [3]; %from 2-16
					case ['7M1']
						pz2go = [16]; %from 10-16
						m12go = [1]; %from 1-10
					case ['7Q1']
						pz2go = [16];%from 2-4, 6-11, 13-16
						m12go = [7]; %from 7-16
					otherwise
						error(['no rat data of ' thisRATID])
				end %thisRATID
			elseif ageRUN == 2
				switch thisRATID
					case ['1B3']
						pz2go = [14];
						m12go = [3];
					case ['6N1']
						pz2go = [14];
						m12go = [1];
					case ['7H2']
						pz2go = [16];
						m12go = [4];
					case ['7M2']
						pz2go = [10];
						m12go = [1];
					case ['7N1']
						pz2go = [15];
						m12go = [1];
					case ['7P1']
						pz2go = [16];
						m12go = [1];
					case ['7Q1']
						pz2go = [16];
						m12go = [2];
					case ['8B2']
						pz2go = [9];
						m12go = [1];
					case ['8K1']
						pz2go = [7];
						m12go = [8];
					otherwise
						error(['no rat data of ' thisRATID])
				end %thisRATID
			end %ageRUN == 1
			targSaveROOT = [preprocROOT '\' ageGROUP{ageRUN} '\' thisRATID]; if ~exist(targSaveROOT) mkdir(targSaveROOT); end

			load([datROOT '\' ageGROUP{ageRUN} '\' ratDat(ratRUN).name]); 	%loading raw data
				if strcmp(thisRATID, '7H2') == 1
					clear Def_*;
					load(['Y:\Jin\Data\middh\P12_7H2_NEWRANGES.mat']);
				end %strcmp(thisRATID, '7H2') == 1
			if ~exist([targSaveROOT '\preproc.mat'])
				preproc = [];
			else
				load([targSaveROOT '\preproc.mat'])
			end %~exist([targSaveROOT '\preproc.mat'])

			%Sleep states time information extract/parsing
			if isfield(preproc, 'slTime') == 0
				slTime = [];
				var2eval = whos(['*Def_*']);
				for sleepRUN = 1:1:size(sleepStage, 1) %AS, QS, then WA
					if contains(var2eval(sleepRUN).name, sleepStage{sleepRUN})
						thisSleepStage = eval([var2eval(sleepRUN).name]);

						thisSLTime = [];
						thisSLTime = [thisSleepStage.times(1:2:(size(thisSleepStage.times, 1) - 1)) thisSleepStage.times(2:2:(size(thisSleepStage.times, 1)))];
						slTime{sleepRUN, 1} = thisSLTime;
					else
						error(['check the raw file']);
					end %contains(var2eval(sleepRUN).name, sleepStage{sleepRUN})
				end %sleepRUN = 1:1:size(sleepStage, 1)

				preproc.slTime = slTime;
				save([targSaveROOT '\preproc.mat'], 'preproc', '-v7.3'); clear var2eval Def_*;
			else
				slTime = preproc.slTime;
			end %isfield(preproc, 'slTime') == 0

			%respiratory signal & respiratory time space
			if isfield(preproc, 'sResp') == 0
				sResp = [];
				tSpace = [];

				try; sResp = Respiration; catch; sResp = resp_1; end
				tSpace = sResp.times;
				sResp = sResp.values;

				preproc.tSpace = tSpace;
				preproc.sResp = sResp;
				save([targSaveROOT '\preproc.mat'], 'preproc', '-v7.3');
			else
				sResp = preproc.sResp;
				tSpace = preproc.tSpace;
			end %isfield(preproc, 'sResp') == 0

			%lfp signals & lfp time space
			if isfield(preproc, 'pzlfp') == 0
				pzlfp = [];
				m1lfp = [];
				lfpTSpace = [];
					lfpTSpacePZ = [];
					lfpTSpaceM1 = [];

				%PZ
				pzVar2eval = whos(['PZ_LFP*']);
				for lfpRUN = 1:1:size(pzVar2eval, 1)
					tempName = pzVar2eval(lfpRUN).name;
					whereUnderscore = find(tempName == '_');
					if tempName(whereUnderscore(end)+1:end) == num2str(pz2go)
						lfp2use = eval([pzVar2eval(lfpRUN).name]);

						pzlfp = lfp2use.values;
						lfpTSpacePZ = lfp2use.times;
						break;
					end %tempName(whereUnderscore(end)+1:end) == num2str(pz2go)
					clear tempName;
				end %lfpRUN = 1:1:size(pzVar2eval, 1)

				pzTIdx = dsearchn(lfpTSpacePZ, [tSpace(1) ; tSpace(end)]);
				pzlfp = pzlfp(pzTIdx(1):pzTIdx(2)); clear lfpTSpacePZ pzTIdx;

				%M1
				m1Var2eval = whos(['Cortical*LFP*']);
				for lfpRUN = 1:1:size(m1Var2eval, 1)
					tempName = m1Var2eval(lfpRUN).name;
					whereUnderscore = find(tempName == '_');
					if tempName(whereUnderscore(end)+1:end) == num2str(m12go)
						lfp2use = eval([m1Var2eval(lfpRUN).name]);

						m1lfp = lfp2use.values;
						lfpTSpaceM1 = lfp2use.times;
						break;
					end %tempName(whereUnderscore(end)+1:end) == num2str(m12go)
					clear tempName;
				end %lfpRUN = 1:1:size(m1Var2eval, 1)

				m1TIdx = dsearchn(lfpTSpaceM1, [tSpace(1) ; tSpace(end)]);
				m1lfp = m1lfp(m1TIdx(1):m1TIdx(2)); clear lfpTSpaceM1 m1TIdx;

				if length(pzlfp) ~= length(m1lfp)
					error(['check raw lfp signals']);
				end %length(pzlfp) ~= length(m1lfp)
				lfpTSpace = linspace(tSpace(1), tSpace(end), length(pzlfp));

				preproc.pzlfp = pzlfp;
				preproc.m1lfp = m1lfp;
				preproc.lfpTSpace = lfpTSpace;
				save([targSaveROOT '\preproc.mat'], 'preproc', '-v7.3'); clear pzVar2eval m1Var2eval PZ_LFP* Cortical_LFP*;
			else
				pzlfp = preproc.pzlfp;
				m1lfp = preproc.m1lfp;
				lfpTSpace = preproc.lfpTSpace;
			end %isfield(preproc, 'pzlfp') == 0

			%extract midha time
			if isfield(preproc, 'midhaWA') == 0
				midhaAS = zeros(1, length(tSpace)); midhaQS = zeros(1, length(tSpace)); midhaWA = zeros(1, length(tSpace));
				midhaASraw = zeros(1, length(lfpTSpace)); midhaQSraw = zeros(1, length(lfpTSpace)); midhaWAraw = zeros(1, length(lfpTSpace));

				for sleepRUN = 1:1:size(sleepStage, 1) 	%AS, QS, then WA
					thisTS = slTime{sleepRUN}; thisState2use = zeros(1, length(tSpace)); thisState2useRaw = zeros(1, length(lfpTSpace));

					for tsRUN = 1:1:size(thisTS, 1)
						tsIdx2go = dsearchn(tSpace, transpose(thisTS(tsRUN, :)));
						thisState2use(1, tsIdx2go(1):tsIdx2go(2)) = 1; clear tsIdx2go;

						tsIdx2go = dsearchn(transpose(lfpTSpace), transpose(thisTS(tsRUN, :)));
						thisState2useRaw(1, tsIdx2go(1):tsIdx2go(2)) = 1; clear tsIdx2go;
					end %tsRUN = 1:1:size(thisTS, 1)

					if sleepRUN == 1
						midhaAS = thisState2use;
						midhaASraw = thisState2useRaw;
					elseif sleepRUN == 2
						midhaQS = thisState2use;
						midhaQSraw = thisState2useRaw;
					elseif sleepRUN == 3
						midhaWA = thisState2use;
						midhaWAraw = thisState2useRaw;
					end %sleepRUN == 1

					clear thisTS;
				end %sleepRUN = 1:1:size(sleepStage, 1)
				%figure; subplot(1, 2, 1); imagesc([midhaAS ; midhaQS ; midhaWA]); subplot(1, 2, 2); imagesc([midhaASraw ; midhaQSraw ; midhaWAraw]);

				preproc.midhaAS = midhaAS;
				preproc.midhaQS = midhaQS;
				preproc.midhaWA = midhaWA;
				preproc.midhaASraw = midhaASraw;
				preproc.midhaQSraw = midhaQSraw;
				preproc.midhaWAraw = midhaWAraw;
				save([targSaveROOT '\preproc.mat'], 'preproc', '-v7.3');
			else
				midhaAS = preproc.midhaAS;
				midhaQS = preproc.midhaQS;
				midhaWA = preproc.midhaWA;
				midhaASraw = preproc.midhaASraw;
				midhaQSraw = preproc.midhaQSraw;
				midhaWAraw = preproc.midhaWAraw;
			end %isfield(preproc, 'midhaWA') == 0

			%get IBR
			if isfield(preproc, 'islerQS') == 0
				%IBR & Isler QS computions were done with WA-truncated respiratory signals
				thisSResp = sResp;
				thisSResp(midhaWA == 1) = nan;
				[thisPks thisPksLocs] = findpeaks(thisSResp, tSpace, 'MinPeakDistance', .25, 'MinPeakProminence', .00001); %thisPksLocs in s resolution

				b2bInterval = diff(thisPksLocs);
				thisIBR = ones(length(b2bInterval), 1) ./ b2bInterval; thisIBR(thisIBR < .5) = nan;
				thisIBR = [nan ; thisIBR];

				ibrTSpace = [min(tSpace):1:max(tSpace)]; varMat = nan(1, length(ibrTSpace));
				for tRUN = 1:1:length(ibrTSpace)
					tIdx2use = dsearchn(thisPksLocs, ibrTSpace(tRUN));
					if abs(thisPksLocs(tIdx2use) - ibrTSpace(tRUN)) < 1
						minTIdx = dsearchn(thisPksLocs, (ibrTSpace(tRUN) - 1));
						maxTIdx = dsearchn(thisPksLocs, (ibrTSpace(tRUN) + 1));

						period2use = [thisIBR(minTIdx:maxTIdx)];
						varMat(tRUN) = nanvar(period2use);
					end %abs(thisPksLocs(tIdx2use) - ibrTSpace(tRUN)) < 1
					clear tIdx2use minTIdx maxTIdx period2use;
				end %tRUN = 1:1:length(ibrTSpace)
				zVarMat = (varMat - nanmean(varMat, 2)) ./ nanstd(varMat); %hist(zVarMat, 50); <- check the positve skewness
					%zVarMat = fillmissing(zVarMat, 'movmedian', round(SFResp));
	            
	            islerQS = zeros(1, length(tSpace)); islerAS = zeros(1, length(tSpace));
	            for tRUN = 1:1:length(tSpace)
	                thisTTT = dsearchn(transpose(ibrTSpace), tSpace(tRUN));
	                if abs(ibrTSpace(thisTTT) - tSpace(tRUN)) < .5
	                    if zVarMat(thisTTT) < prctile(zVarMat, 75)
	                        islerQS(tRUN) = 1; islerAS(tRUN) = 0;
	                    else
	                        islerQS(tRUN) = 0; islerAS(tRUN) = 1;
	                    end %zVarMat(thisTTT) < prctile(zVarMat, 75)
	                end %abs(ibrTSpace(thisTTT) - tSpace(tRUN)) < 1
	            end %tRUN = 1:1:length(tSpace)

	            islerQS(midhaWA == 1) = 0; islerAS(midhaWA == 1) = 0;

	            preproc.islerQS = islerQS;
	            preproc.islerAS = islerAS;
	            save([targSaveROOT '\preproc.mat'], 'preproc', '-v7.3');
			else
				islerQS = preproc.islerQS;
				islerAS = preproc.islerAS;
			end %isfield(preproc, 'islerQS') == 0

			if isfield(preproc, 'islerQSraw') == 0
				islerQSraw = zeros(1, length(lfpTSpace));

				go = true; idx2dig = 1;
				while go
					if islerQS(idx2dig) == 0
						idx2dig = idx2dig + 1;
					elseif islerQS(idx2dig) == 1
						where2start = idx2dig;
						where2end = where2start + 1;
						go2 = true;
						while go2
							if islerQS(where2end) == 1
								where2end = where2end + 1;
	                        elseif islerQS(where2end) == 0
	                            where2end = where2end - 1; go2 = false;
	                        elseif isnan(islerQS)
	                        end %islerQS(where2end) == 1
	                        
	                        if (where2end >= length(islerQS)) | (where2start >= length(islerQS))
	                            go2 = false;
	                        end %(where2end >= length(islerQS)) | (where2start >= length(islerQS))
	                    end %go2
	                    rawTSIdx = dsearchn(transpose(lfpTSpace), [tSpace(where2start) ; tSpace(where2end)]);
	                    islerQSraw(rawTSIdx(1):rawTSIdx(2)) = 1;

	                    idx2dig = where2end + 1; clear where2start where2end;
					end %islerQS(idx2dig) == 0
					if idx2dig >= length(islerQS)
						go = false;
					end %idx2dig >= length(islerQS)
				end %go
				preproc.islerQSraw = islerQSraw;
				save([targSaveROOT '\preproc.mat'], 'preproc', '-v7.3');
			else
				islerQSraw = preproc.islerQSraw;
			end %isfield(preproc, 'islerQSraw') == 0

			if isfield(preproc, 'islerASraw') == 0
				islerASraw = zeros(1, length(lfpTSpace));

				go = true; idx2dig = 1;
				while go
					if islerAS(idx2dig) == 0
						idx2dig = idx2dig + 1;
					elseif islerAS(idx2dig) == 1
						where2start = idx2dig;
						where2end = where2start + 1;
						go2 = true;
						while go2
							if islerAS(where2end) == 1
								where2end = where2end + 1;
	                        elseif islerAS(where2end) == 0
	                            where2end = where2end - 1; go2 = false;
	                        elseif isnan(islerAS)
	                        end %islerAS(where2end) == 1
	                        
	                        if (where2end >= length(islerAS)) | (where2start >= length(islerAS))
	                            go2 = false;
	                        end %(where2end >= length(islerAS)) | (where2start >= length(islerAS))
	                    end %go2
	                    rawTSIdx = dsearchn(transpose(lfpTSpace), [tSpace(where2start) ; tSpace(where2end)]);
	                    islerASraw(rawTSIdx(1):rawTSIdx(2)) = 1;

	                    idx2dig = where2end + 1; clear where2start where2end;
					end %islerAS(idx2dig) == 0
					if idx2dig >= length(islerAS)
						go = false;
					end %idx2dig >= length(islerAS)
				end %go
				preproc.islerASraw = islerASraw;
				save([targSaveROOT '\preproc.mat'], 'preproc', '-v7.3');
			else
				islerASraw = preproc.islerASraw;
			end %isfield(preproc, 'islerASraw') == 0

			if isfield(preproc, 'resp2SDraw') == 0
				sdCri = 2;
				pwResp = abs(sResp) .^ 2;
					meanResp = nanmean(pwResp, 1); sdResp = nanstd(pwResp, [], 1);

				resp2SD = ones(1, length(pwResp)); resp2SDraw = ones(1, length(lfpTSpace));
				outlierIdx = find((pwResp > (meanResp + sdResp * sdCri)) | (pwResp < (meanResp - sdResp * sdCri)) == 1);

				for outlierRUN = 1:1:length(outlierIdx)
					preIdx = outlierIdx(outlierRUN) - 1; %SF = 30, so 1 bin > 33.3333ms
					postIdx = outlierIdx(outlierRUN) + 1;

					resp2SD(preIdx:postIdx) = 0;

					%under high SF
					remNBin = 20; %20ms before and after the outlier points
	                highSFidx2use = dsearchn(transpose(lfpTSpace), [tSpace(outlierIdx(outlierRUN))]);

	                resp2SDraw((highSFidx2use(1) - remNBin):(highSFidx2use(1) + remNBin)) = 0;
				end %outlierRUN = 1:1:length(outlierIdx)

				preproc.resp2SDraw = resp2SDraw;
	            preproc.resp2SD = resp2SD;
	            save([targSaveROOT '\preproc.mat'], 'preproc', '-v7.3');
			else
				resp2SDraw = preproc.resp2SDraw;
	            resp2SD = preproc.resp2SD;
			end %isfield(preproc, 'resp2SDraw') == 0

			%get raw delta lfps
			if isfield(preproc, 'pzDelta') == 0
				pzDelta = [];
				m1Delta = [];

				pzDelta = jkBandpassFilt(pzlfp, SFtarg, deltaFreq);
				m1Delta = jkBandpassFilt(m1lfp, SFtarg, deltaFreq);

				preproc.pzDelta = pzDelta;
				preproc.m1Delta = m1Delta;
				save([targSaveROOT '\preproc.mat'], 'preproc', '-v7.3');
			else
				%pzDelta = preproc.pzDelta;
				%m1Delta = preproc.m1Delta;
			end %isfield(preproc, 'pzDelta') == 0

			if isfield(preproc, 'pzDelta1MedRaw') == 0
				[pzDelta1MedRaw ts2go] = jkParsePowerSection(pzDelta, lfpTSpace, SFtarg, 1, 3, 1); %1 med, 3 sec, + smoothing!
				pzDelta1Med = zeros(1, length(tSpace));
				for tRUN = 1:1:size(ts2go, 1)
					tIdx2use = dsearchn(tSpace, transpose(ts2go(tRUN, :)));
					pzDelta1Med(tIdx2use(1):tIdx2use(2)) = 1;
					clear tIdx2use;
				end %tRUN = 1:1:size(ts2go, 1)

				preproc.pzDelta1MedRaw = transpose(pzDelta1MedRaw);
				preproc.pzDelta1Med = pzDelta1Med;
				save([targSaveROOT '\preproc.mat'], 'preproc', '-v7.3');
			else
				pzDelta1MedRaw = preproc.pzDelta1MedRaw;
				pzDelta1Med = preproc.pzDelta1Med;
			end %isfield(preproc, 'pzDelta1MedRaw')

			if isfield(preproc, 'm1Delta1MedRaw') == 0
				[m1Delta1MedRaw ts2go] = jkParsePowerSection(m1Delta, lfpTSpace, SFtarg, 1, 3, 1); %1 med, 3 sec, + smoothing!
				m1Delta1Med = zeros(1, length(tSpace));
				for tRUN = 1:1:size(ts2go, 1)
					tIdx2use = dsearchn(tSpace, transpose(ts2go(tRUN, :)));
					m1Delta1Med(tIdx2use(1):tIdx2use(2)) = 1;
					clear tIdx2use;
				end %tRUN = 1:1:size(ts2go, 1)

				preproc.m1Delta1MedRaw = transpose(m1Delta1MedRaw);
				preproc.m1Delta1Med = m1Delta1Med;
				save([targSaveROOT '\preproc.mat'], 'preproc', '-v7.3');
			else
				m1Delta1MedRaw = preproc.m1Delta1MedRaw;
				m1Delta1Med = preproc.m1Delta1Med;
			end %isfield(preproc, 'm1Delta1MedRaw')

			%get raw spike timestamps
			if isfield(preproc, 'pzMUA') == 0
				pzMUA = [];
				pzSC = [];

				mua2eval = whos(['Raw*mua*']);
				for muaRUN = 1:1:size(mua2eval, 1)
					thisMUA = eval([mua2eval(muaRUN).name]);
					pzMUA{muaRUN, 1} = thisMUA.times;

					clear thisMUA;
				end %muaRUN = 1:1:size(mua2eval, 1)

				sc2eval = whos(['Raw*good*']);
				for scRUN = 1:1:size(sc2eval, 1)
					thisSC = eval([sc2eval(scRUN).name]);
					pzSC{scRUN, 1} = thisSC.times;

					clear thisSC;
				end %scRUN = 1:1:size(sc2eval, 1)

				preproc.pzMUA = pzMUA;
				preproc.pzSC = pzSC;
				save([targSaveROOT '\preproc.mat'], 'preproc', '-v7.3'); clear mua2eval sc2eval Raw*mua* Raw*good*;
			else
				%pzMUA = preproc.pzMUA;
				%pzSC = preproc.pzSC;
			end %isfield(preproc, 'pzMUA') == 0

			%get rasters under 1ms
			if isfield(preproc, 'pzRasters') == 0
				tic;
				pzRasters = []; cellCount = 1;
	%				allSpkTS = [];
	%			pzIFR = []; 	%smoohted data of sum of pzRasters

				for ctypeRUN = 1:1:2 	%mua then sc
					if ctypeRUN == 1
						thisTS2play = pzMUA;
					else
						thisTS2play = pzSC;
					end %ctypeRUN == 1

					for cell2RUN = 1:1:size(thisTS2play, 1)
						thisRaster = zeros(1, length(lfpTSpace));
						thisRaster(dsearchn(transpose(lfpTSpace), thisTS2play{cell2RUN})) = 1;

						pzRasters = [pzRasters ; thisRaster];
						cellCount = cellCount + 1;
					end %cell2RUN = 1:1:size(thisTS2play, 1)
					clear thisTS2play;
				end %ctypeRUN = 1:1:2 	%mua then sc

				preproc.pzRasters = pzRasters;
				preproc.pzIFR = nansum(pzRasters, 1);
				save([targSaveROOT '\preproc.mat'], 'preproc', '-v7.3');
				toc;
			else
				%pzRasters = preproc.pzRasters;
				%pzIFR = nansum(pzRasters, 1);
			end %isfield(preproc, 'pzRasters') == 0

			%version#1
			nFIGROW = 2; nFIGCOL = 3; txtSF = {['low'] ; ['high']};
				txtXTick2go{1, 1} = {['Resp2SD'], ['MidhaAS'], ['IslerAS']};
				txtXTick2go{2, 1} = {['Resp2SD'], ['MidhaQS'], ['IslerQS'], ['PZ Delta 1Med'], ['M1 Delta 1Med']};
				txtXTick2go{3, 1} = {['Resp2SD'], ['MidhaWA']};
			picID = figure('Color', 'w', 'Position', [25 25 1200 680]);
			for sfRUN = 1:1:2
				if sfRUN == 1
					nDats = length(tSpace);

					thisR2go = resp2SD;
				elseif sfRUN == 2
					nDats = length(lfpTSpace);

					thisR2go = resp2SDraw;
				end %sfRUN == 1

				for slRUN = 1:1:3
					if slRUN == 1
						if sfRUN == 1
							thisStates = midhaAS;
							thisIsler = islerAS;
						elseif sfRUN == 2
							thisStates = midhaASraw;
							thisIsler = islerASraw;
						end %sfRUN == 1

						remDatSpace = [thisR2go ; thisStates ; thisIsler];
					elseif slRUN == 2
						if sfRUN == 1
							thisStates = midhaQS;
							thisIsler = islerQS;
							thisPZMed = pzDelta1Med;
							thisM1Med = m1Delta1Med;
						elseif sfRUN == 2
							thisStates = midhaQSraw;
							thisIsler = islerQSraw;
							thisPZMed = pzDelta1MedRaw;
							thisM1Med = m1Delta1MedRaw;
						end %sfRUN == 1

						remDatSpace = [thisR2go ; thisStates ; thisIsler ; thisPZMed ; thisM1Med];
					elseif slRUN == 3
						if sfRUN == 1
							thisStates = midhaWA;
						elseif sfRUN == 2
							thisStates = midhaWAraw;
						end %sfRUN == 1

						remDatSpace = [thisR2go ; thisStates];
					end %slRUN == 1

					subplot(nFIGROW, nFIGCOL, slRUN + (sfRUN - 1) * nFIGCOL, 'FontSize', szTXT);
					hold on;
						for rowRUN = 1:1:size(remDatSpace, 1)
							plot(rowRUN, nansum(nansum(remDatSpace(1:rowRUN, :), 1) == rowRUN) / nDats * 100, 'ks');
						end %rowRUN = 1:1:size(remDatSpace, 1)
						plot([0 100], ones(1, 2) .* 50, ':k');
						plot([0 100], ones(1, 2) .* 10, ':k');
						plot([0 100], ones(1, 2) .* 5, ':r');
					hold off;
					set(gca, 'XLim', [.5 rowRUN + .5], 'XTick', 1:rowRUN, 'XTickLabel', txtXTick2go{slRUN}, 'YLim', [0 100], 'YTick', 0:25:100);
					if slRUN == 1 ylabel(['prop of available data']); end
					title({[thisRATID '-' txtSF{sfRUN} ' SF'] ; [sleepStage{slRUN}]});

					clear remDatSpace;
				end %slRUN = 1:1:3
			end %sfRUN = 1:1:2
			saveas(picID, [visSum '\R' thisRATID '-SanityCheck.bmp']); close all;

			%version#2
			nFIGROW = 1; nFIGCOL = 2; slRUN = 2; %QS only
			txt1Min = {['No1MinCri']; ['Yes1MinCri']}; txt1MinId = 1;
			thisPropVal2put = [];
			picID = figure('Color', 'w', 'Position', [25 25 600 480]);
			for sfRUN = 1:1:2
				if sfRUN == 1
					thisR2go = resp2SD;
					thisStates = midhaQS;
					thisIsler = islerQS;
					thisPZMed = pzDelta1Med;
					thisM1Med = m1Delta1Med;

					thisSF = 30;
				elseif sfRUN == 2
					thisR2go = resp2SDraw;
					thisStates = midhaQSraw;
					thisIsler = islerQSraw;
					thisPZMed = pzDelta1MedRaw;
					thisM1Med = m1Delta1MedRaw;

					thisSF = 1000;
				end %sfRUN == 1

				remDatSpace = [thisR2go ; thisStates ; thisIsler ; thisPZMed ; thisM1Med];
				nDats = nansum(nansum(remDatSpace(1:3, :), 1) == 3, 2);
				if nansum(nansum(remDatSpace, 1) == size(remDatSpace, 5), 2) >= (thisSF * 60)
					txt1MinId = 2;
				end %nansum(nansum(remDatSpace, 1) == size(remDatSpace, 5), 2) >= (thisSF * 60)
				thisPropVal2put = [thisPropVal2put txt1MinId];

				subplot(nFIGROW, nFIGCOL, sfRUN, 'FontSize', szTXT);
				hold on;
					for rowRUN = 4:size(remDatSpace, 1)
						%bar(rowRUN - 3, nansum(nansum(remDatSpace(1:rowRUN, :), 1) == rowRUN) / nDats * 100);
						plot(rowRUN - 3, nansum(nansum(remDatSpace(1:rowRUN, :), 1) == rowRUN) / nDats * 100, 'ks');
						thisPropVal2put = [thisPropVal2put nansum(nansum(remDatSpace(1:rowRUN, :), 1) == rowRUN) / nDats * 100];
					end %rowRUN = 4:size(remDatSpace, 1)
					thisYLim = get(gca, 'YLim');
				hold off;
				set(gca, 'XLim', [.5 2.5], 'XTick', 1:2, 'XTickLabel', {['PZ > 1*Med'], ['M1 > 1*Med']}, 'YLim', [min(thisYLim) - 5 max(thisYLim) + 5]);
				if sfRUN == 1 ylabel(['prop to Isler QS (100%)']); end
				title({[thisRATID '-' txtSF{sfRUN} ' SF'] ; [sleepStage{slRUN} ' ' txt1Min{txt1MinId}]});

				clear remDatSpace; txt1MinId = 1;
			end %sfRUN = 1:1:2
			qsPropTab = [qsPropTab ; thisPropVal2put];
			saveas(picID, [visSum '\R' thisRATID '-SanityCheckV2.bmp']); close all;

			clear thisRATID targSaveROOT preproc slTime sResp tSpace pzlfp m1lfp lfpTSpace midhaAS midhaQS midhaWA midhaASraw midhaQSraw midhaWAraw islerQS islerAS islerQSraw islerASraw resp2SDraw resp2SD pzDelta m1Delta pzDelta1MedRaw pzDelta1Med m1Delta1MedRaw m1Delta1Med pzMUA pzSC pzRasters pzIFR;
		end %ratRUN = 1:1:size(ratDat, 1)
		clear ratDat;
	end %ageRUN = 2:1:size(ageGROUP, 1)
	save([visSum '\qsPropTab.mat'], 'qsPropTab');
else
	load([visSum '\qsPropTab.mat']);
end %~exist([visSum '\qsPropTab.mat'])

%version#2 pop
nFIGROW = 1; nFIGCOL = 2; slRUN = 2; szTXT = 7.5; txtSF = {['low'] ; ['high']};%QS only
txt1Min = {['No1MinCri']; ['Yes1MinCri']}; txt1MinId = 1;
picID = figure('Color', 'w', 'Position', [25 25 600 480]);
for sfRUN = 1:1:2
	thisTab2vis = qsPropTab(:, (1:3) + (sfRUN - 1) * 3);

	subplot(nFIGROW, nFIGCOL, sfRUN, 'FontSize', szTXT);
	hold on;
		for indRUN = [1 3 4 5 6 7 8 9]
			if thisTab2vis(indRUN, 1) == 2
				plot([1 2], thisTab2vis(indRUN, 2:3), '.', 'Color', ones(1, 3) .* .75);
			end %thisTab2vis(indRUN, 1) == 2
		end %indRUN = 1:1:size(thisTab2vis, 1)
		errorbar([1 2], nanmean(thisTab2vis(thisTab2vis(:, 1) == 2, 2:3), 1), nanstd(thisTab2vis(thisTab2vis(:, 1) == 2, 2:3), [], 1) ./ sqrt(8), 'ks');
	hold off
	set(gca, 'XLim', [.5 2.5], 'XTick', 1:2, 'XTickLabel', {['PZ > 1*Med'], ['M1 > 1*Med']});
	if sfRUN == 1 ylabel(['prop to Isler QS (100%)']); end
	title({['Pop ' txtSF{sfRUN} ' SF'] ; [sleepStage{slRUN}]});

	clear remDatSpace; txt1MinId = 1;
end %sfRUN = 1:1:2
print -depsc -tiff -r300 -painters 'Pop-Sanity.eps';

[h1 p1 ci1 stats1] = ttest2(qsPropTab([1 3 4 5 6 7 8 9], 2), qsPropTab([1 3 4 5 6 7 8 9], 3))
[h2 p2 ci2 stats2] = ttest2(qsPropTab([1 3 4 5 6 7 8 9], 5), qsPropTab([1 3 4 5 6 7 8 9], 6))