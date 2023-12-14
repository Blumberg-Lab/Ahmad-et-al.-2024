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
	plvROOT = [alzROOT '\XCorrSpkSpkV1']; if ~exist(plvROOT) mkdir(plvROOT); end
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
if ~exist([plvROOT '\xcorrTab.mat'])
	xcorrTab = [];

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

			nFIGROW = 2; nFIGCOL = 2;
            for slRUN = 1:1:2 	%AS & QS
                if slRUN == 1
                    stateDef = [resp2SDraw ; midhaASraw ; islerASraw];
                elseif slRUN == 2
                    stateDef = [resp2SDraw ; midhaQSraw ; islerQSraw ; pzDelta1MedRaw];
                end %slRUN == 1
                dat2go = find(nansum(stateDef, 1) == size(stateDef, 1));

                spkspkDat = pzRasters(:, dat2go);

                %binning with 5ms
                spkspkDatbinned = [];
                szBin = 5; thisQuotient = floor(size(spkspkDat, 2) / szBin);
                for binRUN = 1:1:thisQuotient
                    stBin = 1 + (binRUN - 1) * szBin;
                    edBin = binRUN * szBin;
                    
                    spkspkDatbinned(:, binRUN) = nansum(spkspkDat(:, stBin:edBin), 2);
                end %binRUN = 1:1:thisQuotient

                %chunking with 2s
                spkspkDatbc = []; szChunk = 2000 / szBin;
                for cellRUN = 1:1:size(spkspkDatbinned, 1)
                    thisDat2chk = spkspkDatbinned(cellRUN, :);
                    thisQuotient2 = floor(size(thisDat2chk, 2) / szChunk);
                    thisDat2chk = transpose(thisDat2chk(1, 1:(thisQuotient2 * szChunk)));
                    spkspkDatbc = [spkspkDatbc ; ... 
                                    transpose(nanmean(reshape(thisDat2chk, szChunk, thisQuotient2), 2))];

                    clear thisDat2chk thisQuotient2 thisDat2chk;
                end %cellRUN = 1:1:size(spkspkDatbinned, 1)
                
                corrmap = nan(size(spkspkDatbc, 1), size(spkspkDatbc, 1));
                for ii = 1:1:size(spkspkDat, 1)
                    for jj = 1:1:size(spkspkDat, 1)
                         [dummyR dummyP] = corrcoef(spkspkDatbc(ii, :), spkspkDatbc(jj, :));
                         corrmap(ii, jj) = dummyR(1, 2);
                    end %jj = 1:1:size(spkspkDat, 1)
                end %ii = 1:1:size(spkspkDat, 1)
                xcorrTab{ratRUN, slRUN} = corrmap;

                clear dat2go spkspkDat corrmap;
            end %slRUN = 1:1:2 	%AS & QS

			clear thisRATID preprocLOADROOT preproc slTime sResp tSpace pzlfp m1lfp lfpTSpace midhaAS midhaQS midhaWA midhaASraw midhaQSraw midhaWAraw islerQS islerAS islerQSraw islerASraw resp2SDraw resp2SD pzDelta m1Delta pzDelta1MedRaw pzDelta1Med m1Delta1MedRaw m1Delta1Med pzMUA pzSC pzRasters pzIFR;
		end %ratRUN = 1:1:size(ratIDs2go, 1)
	end %ageRUN = 2:1:size(ageGROUP, 1)
	save([plvROOT '\xcorrTab.mat'], 'xcorrTab');
else
	load([plvROOT '\xcorrTab.mat']);
end %~exist([plvROOT '\plvMasterTab.mat'])

nFIGROW = 2; nFIGCOL = nansum(pzPLVFlag, 1) + 1; ratCount = 0;
picID = figure('Color', [1 1 1], 'Position', [65 65 1000 600]);
avgCoeffTab = [];
for ratRUN = 1:1:size(ratIDs2go, 1)
    if pzPLVFlag(ratRUN) == 1
        ratCount = ratCount + 1;
        for slRUN = 1:1:2
            thisMap = triu(xcorrTab{ratRUN, slRUN}, 1);
            thisMap(thisMap == 0) = nan;

            subplot(nFIGROW, nFIGCOL, ratCount + (slRUN - 1) * nFIGCOL, 'FontSize', szTXT);
            hold on;
                hh = imagesc(thisMap);
            hold off; axis off tight;
            set(hh, 'AlphaData', ~isnan(thisMap));
            if slRUN == 1 title({[ratIDs2go{ratRUN}] ; ['AS (n = ' num2str(size(thisMap, 1)) ')']}); else title(['QS']); end
            
            avgCoeffTab(ratCount, slRUN) = nanmean(atanh(thisMap(:)), 1);
        end %slRUN = 1:1:2
    end %pzPLVFlag(ratRUN) == 1
end %ratRUN = 1:1:size(ratIDs2go, 1)
subplot(nFIGROW, nFIGCOL, [nFIGCOL nFIGCOL * 2], 'FontSize', szTXT);
hold on;
    %errorbar([1 2], nanmean(avgCoeffTab, 1), nanstd(avgCoeffTab, [], 1) ./ sqrt(ratCount));
    boxplot(avgCoeffTab);
    for ratRUN = 1:1:size(avgCoeffTab, 1)
        plot([1 2], avgCoeffTab(ratRUN, :), 'k.');
    end %ratRUN = 1:1:size(avgCoeffTab, 1)
hold off;
set(gca, 'XLim', [.5 2.5], 'XTick', 1:1:2, 'XTickLabel', {['AS'], ['QS']}); ylabel(['r-to-z coeff']);
