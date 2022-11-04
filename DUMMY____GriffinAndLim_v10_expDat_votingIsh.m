% Spectrogram function
addpath( 'C:\Users\Lord Photon\Documents\MATLAB\library_repo\library' )
addpath( '/Users/ben/Documents/MATLAB/library_repo' )
addpath('../phaseRecovery_Data');
outputFigsDir='C:\Users\Lord Photon\Documents\MATLAB\time-frequency analysis\outputFigs_DUMMY_Cross_20221104\';

%% Load data and set time-frequency vectors

plotIter=0;
% % define loop parameters
%%%%%%%%%%%%%%%%%%%%%% ZIG ZAG
% clips_l=[20,40,60];
% nptPerWin_l=[16,32,64];
% bwIRF_l=[30:5:70]*1e9;
% interpAmounts=[1,2,4,8,16];
% %
% clipThresh=60;
% nptPerWin=16;
% bwIRF=45e9; % 50 GHz is optimal
%  interpAmount_t=4 ; % For now, make this a power of 2 (or 1)!!

%
% for icl=2:numel(clips_l)
%     for inptWin=1:numel(nptPerWin_l)
%         for ibwIRF=1:numel(bwIRF_l)
%             for iInter=1:numel(interpAmounts)
%
%
% clipThresh=clips_l(icl);
% nptPerWin=nptPerWin_l(inptWin);
% bwIRF=bwIRF_l(ibwIRF); % 50 GHz is optimal
%  interpAmount_t=interpAmounts(iInter) ; % For now, make this a power of 2 (or 1)!!




%% RTO ZigZag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('RTOZigZag.mat');lowerClip=max(max(spgm))/400;
% load('/Users/ben/Documents/MATLAB/timeFrequencyAnalysis/phaseRecovery_Data/RTOzigzag_deconv.mat');
% load('C:\Users\Lord Photon\Documents\MATLAB\time-frequency analysis\PhaseRecoveryAlgorithms_repo\phaseRecovery_Data/RTOZigZag.mat');
% load('RTOzigzag_deconv.mat');
% load('RTOzigzag_deconv_20221028.mat');
% tIndsExpInterest=2:179;
% fSpecGHz=fSpecGHz;%/(56.4e9*TargetResolution)
% winLen_t=200e-12
% phaseAnalysis=1;
% nFreqElem=numel(fSpecGHz);%2000;
% freqInds=round(linspace(1,numel(fSpecGHz),nFreqElem));
%
%
% plotFilt=0;



%% OSO Cross
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('OsoDeconv_narrow_20221028.mat');spgm=circshift(spgm,-82).^2.5;;%lowerClip=max(max(spgm))/60;
% load('OSO_deconvWithBackground.mat');lowerClip=max(max(spgm))/25;
% load('OSO_deconvNoBackground.mat');spgm=circshift(spgm,-82).^1.5;
% load('OsoDeconv_narrow_irf14p5_40.mat');%.^1.5;
% load('OSOdataCross_filtv2.mat');lowerClip=max(max(spgm))/25;
% load('OSOdataCross.mat');lowerClip=max(max(spgm))/25;

clips_l=[20,70,140];
fns={'OsoDeconv_narrow_irf11p5_26','OsoDeconv_narrow_irf11p5_40','OsoDeconv_narrow_irf13_50',...
    'OsoDeconv_narrow_irf14p5_40'}
nptPerWin_l=[64,128];
bwIRF_l=[400,  800,  8000]*1e9;
interpAmounts=[1,4,8];



for icl=1:numel(clips_l)
    for ifn=1:numel(fns);
        for inptWin=1:numel(nptPerWin_l)
            for ibwIRF=1:numel(bwIRF_l)
                for iInter=1:numel(interpAmounts)
                    itStr=num2str([ icl, ifn, inptWin, ibwIRF, iInter])
                    
                    load(fns{ifn});
                    clipThresh=clips_l(icl);
                    nptPerWin=nptPerWin_l(inptWin);
                    bwIRF=bwIRF_l(ibwIRF); % 50 GHz is optimal
                    interpAmount_t=interpAmounts(iInter) ; % For now, make this a power of 2 (or 1)!!
                    
                    
                    
                    tIndsExpInterest=65:210;
                    winLen_t=62.5e-12
                    
                    phaseAnalysis=2;
                    fMaxStft=500e9;
                    mIRF=8
                    
                    % clipThresh=80;
                    % nptPerWin=64;
                    % bwIRF=800e9; % 50 GHz is optimal
                    % interpAmount_t=2 ; % For now, make this a power of 2 (or 1)!!
                    
                    
                    
                    
                    
                    %% voting & other algorithm parameters
                    maxIteration=30; % Max number of iteration in GLA
                    % interpAmount_t=4 ; % For now, make this a power of 2 (or 1)!!
                    interpAmount_f=1  ; % For now, make this a power of 2 (or 1)!!
                    
                    % algorithm filter parameters
                    fMaxStftFilter=fSpecGHz(end)*1e9*2.5;
                    filtm=10;
                    plotFilt=0;
                    
                    % voting parameters
                    nIter=1;
                    maxIterationVoted=40;
                    
                    
                    
                    %% Filter signal (i.e., reconvolve)
                    
                    dtSpgmRaw=winLen_t/numel(spgm(:,1));
                    tRaw=(1:numel(dtSpgmRaw))*dtSpgmRaw;
                    fRaw=linspace(-1/(2*dtSpgmRaw),1/(2*dtSpgmRaw),numel(spgm));
                    yReshaped=reshape(spgm,[1,numel(spgm)]);
                    yReconv=abs(filterSig(mIRF,bwIRF,yReshaped,1,fRaw));
                    spgm=abs(reshape(yReconv,size(spgm)));
                    
                    
                    spgmIni=spgm;%.^0.8;%(freqInds,:);
                    
                    %% restrict time axis and clip lower values to avoid noise issues.
                    spgmExpRaw=spgmIni(:,tIndsExpInterest);
                    lowerClip=max(max(spgm))/clipThresh;
                    spgmExpRaw(spgmExpRaw<lowerClip)=0;%
                    spgmExpRaw(spgmExpRaw>lowerClip)=spgmExpRaw(spgmExpRaw>lowerClip)-lowerClip;
                    
                    
                    
                    
                    %% Time frequency vectors definition & interpolation to have integer number of lenses
                    
                    fSpecExp=fSpecGHz*1e9;  tSpecExp=tSpecns(tIndsExpInterest)*1e-9;
                    tWind=numel(tSpecExp)*winLen_t;         % The total length of t dictates the minimum required frequency resolution
                    TargetResolution=winLen_t/nptPerWin;    % Target temporal resolution of the STFT (i.e., frequency span)
                    
                    lent=round(tWind/TargetResolution);
                    
                    dt=tWind/lent;tWind=lent*dt;
                    t=(1:lent)*dt;
                    Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
                    fG=f*10^-9;tps=t*10^12;
                    
                    
                    fspgm_raw=f;
                    tspgm_raw=tSpecExp;
                    SUT=zeros(1,lent);
                    
                    % Zero pad if needed
                    if 1/(2*TargetResolution)>fSpecExp(end)
                        
                        [a,nZeros]=find(f>-abs(fSpecGHz(end)*1e9),1); nZeros=nZeros-1;
                        f_tls1=[f(1:nZeros-1),fSpecExp,-flip(f(1:nZeros-1))];
                        padded1=[zeros(nZeros-1,numel(spgmExpRaw(1,:)));spgmExpRaw;zeros(nZeros-1,numel(spgmExpRaw(1,:)))];
                        spgmExpRawPadded=interp2fun(tspgm_raw,f_tls1,padded1,tspgm_raw,f);
                        spgmExpRawPadded(1:nZeros,:)=0;spgmExpRawPadded(lent-nZeros:end,:)=0;
                        
                        % Canot use imresize since interpolation querie points are
                        % irregularly spaced
                        %     spgmExpRawPadded=imresize(padded1,[numel(f), numel(tspgm_raw)]);
                        %     spgmExpRawPadded(1:nZeros,:)=0;spgmExpRawPadded(lent-nZeros:end,:)=0;
                    else
                        spgmExpRawPadded=spgmExpRaw;
                    end
                    
                    
                    
                    
                    %% window Setup
                    
                    winLen=winLen_t/dt;
                    if abs(winLen-round(winLen))<0.01
                        warning(['small rounding of winLen, from ' num2str(winLen_t) ' to ' num2str(winLen*dt)])
                        winLen=round(winLen);
                    end
                    
                    
                    winInc=winLen; % without interpolation, winInc is always equal to winLen
                    winSec=ones(1,winLen);
                    % windowInds=(1:lent)-lent/2;win=superGauss(0,winLen/2,8,windowInds,0);
                    win=zeros(1,lent); win(round(lent/2)-floor(winLen/2):round(lent/2)+ceil(winLen/2)-1)=winSec; % make window as a pure square
                    win=circshift(win,lent/2);
                    
                    nIncs=lent/winInc; % By making winInc a power of 2, we can make sure to have an integer number of windows.
                    windowCenters=(1:nIncs)*winInc;
                    
                    
                    
                    
                    
                    %% Setup interpolation
                    nIncsInterp=nIncs*interpAmount_t;
                    windowCentersInterp=(1:nIncsInterp)*winInc/interpAmount_t;
                    winLenInterp=numel(fspgm_raw)*interpAmount_f;
                    winInterp=interp1(linspace(0,1,lent),win,linspace(0,1,winLenInterp));
                    overlapAmount=interpAmount_t*(sum(win))/winInc;%numel(winIndsInterp)/(windowCentersInterp(2)-windowCentersInterp(1)); % This is the "overlapamount" AFTER interpolation
                    analysisWin=winInterp/overlapAmount; % Analysis window for the inverse spgm
                    
                    tspgm=linspace(tspgm_raw(1),tspgm_raw(end),numel(tspgm_raw)*interpAmount_t);
                    fspgm=linspace(fspgm_raw(1),fspgm_raw(end),numel(fspgm_raw)*interpAmount_f);%fspgm_raw;
                    
                    if interpAmount_t>1
                        spgm=abs(imresize(spgmExpRawPadded,[lent,numel(tspgm)],'bicubic'));%'lanczos3'));
                    else
                        spgm=spgmExpRawPadded;
                    end
                    %
                    % % Remove the extrapolated noise
                    % numToDecimate=round((1/(2*dt)-fSpecExp(end))/df);
                    % spgm(1:numToDecimate,:)=0; spgm(end-numToDecimate:end,:)=0;
                    
                    
                    
                    %% plot of all different spgms
                    
                    
                    %
                    % figure;
                    % subplot(3,1,1)
                    % imagesc(tSpecns,fSpecGHz,spgmIni);
                    % subplot(3,1,2)
                    % imagesc(spgmExpRawPadded)
                    % subplot(3,1,3)
                    % imagesc(tspgm,fspgm,spgm)
                    %
                    %
                    %
                    
                    clear padded1 spgmExpRaw spgmExpRawPadded
                    
                    
                    %% Iterative Griffin and Lim algorithm
                    
                    % initialize output matrices for  voting method
                    SiOut=zeros(numel(spgm(:,1)),numel(spgm(1,:)),nIter);
                    xtOut=zeros(numel(SUT),nIter);
                    
                    for iIter=1:nIter
                        
                        S0=sqrt(spgm).*exp(1j*rand(size(spgm))*2*pi);
                        
                        iniPhase=exp(1j*normrnd(0,pi,size(S0)));
                        
                        xt0=get_istft_fullSigLen(lent,windowCentersInterp,analysisWin,Fs,nIncsInterp,S0.*iniPhase);
                        
                        
                        % Convergence criterion
                        di=zeros(1,maxIteration);
                        diC=di;
                        diR=di;
                        
                        
                        
                        %     if iIter==1; plotIter=1;end
                        
                        [xt,Si,diCPlot,diRPlot]=DUMMY__phaseRecovLoop(nIncsInterp,windowCentersInterp,lent,winInterp,winLen,t,dt,xt0,tspgm,fspgm,spgm,SUT,analysisWin,Fs,maxIteration,{'Inconsistency','frogE'},plotIter,filtm,fMaxStftFilter,plotFilt,f);
                        
                        
                        
                        
                        SiOut(:,:,iIter)=Si;
                        xtOut(:,iIter)=xt;
                        iIter
                    end
                    
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%% FOR PARAMETRIC ANALYSIS
                    alldiCs( icl, ifn, inptWin, ibwIRF, iInter)=diCPlot(end);
                    alldiRs( icl, ifn, inptWin, ibwIRF, iInter)=diRPlot(end);
                    
                    
                    timeNowStr=getTimeString();
                    %                     save_figpng([outputFigsDir 'spgm_' itStr '_' timeNowStr]);
                    f1=gcf;
                    saveas(f1,[outputFigsDir 'spgm_' itStr '_' timeNowStr '.png'],'png')
                    
                    
                    SiOut_Raw=SiOut;
                    
                    
                    
                    if nIter~=1
                        
                        
                        
                        
                        
                        ss=size(SiOut); s1=ss(1); s2=ss(2);
                        for j=1:nIter;
                            SiOut(:,:,nIter)=SiOut(:,:,nIter).*exp(-1j*angle(SiOut(round(s1/2),round(s2/2),nIter))).*exp(1j*angle(angle(SiOut(round(s1/2),round(s2/2),1))));
                        end
                        angle(SiOut(round(s1/2),round(s2/2),:))
                        
                        pairs=nchoosek(1:nIter,2);
                        SiDiffs=zeros(numel(spgm(:,1)),numel(spgm(1,:)),numel(pairs(:,1)));
                        
                        
                        for iVot=1:numel(pairs(:,1))
                            
                            SiDiffs(:,:,iVot)=abs(SiOut(:,:,pairs(iVot,1))-SiOut(:,:,pairs(iVot,2)));
                            
                        end
                        
                        [~,minPair]=min(SiDiffs,[],3);
                        pairs1=pairs(:,1); pairs2=pairs(:,2);
                        allPtsPairs1=pairs1(minPair);allPtsPairs2=pairs2(minPair);
                        
                        % SiMean=mean(SiPairOpt,3);
                        
                        [m,n,k] = size(SiOut);
                        [q,w] = ndgrid(1:m,1:n);
                        Sipair1 = SiOut(sub2ind([m,n,k],q,w,allPtsPairs1));
                        Sipair2 = SiOut(sub2ind([m,n,k],q,w,allPtsPairs2));
                        SiPairOpt=mean(cat(3,Sipair1,Sipair2),3);
                        
                        % SiMean=mean(SiOut,3);
                        xt0=get_istft_fullSigLen(lent,windowCentersInterp,analysisWin,Fs,nIncsInterp,SiPairOpt);
                        plotIter=1;
                        % xt0=mean(xtOut,2).';
                        
                        % xt0=get_istft_fullSigLen(lent,windowCentersInterp,analysisWin,Fs,nIncsInterp,sqrt(spgm).*exp(1j*angle(SiMean)));
                        [xt,Si,diCPlot,diRPlot]=DUMMY__phaseRecovLoop(nIncsInterp,windowCentersInterp,lent,winInterp,winLen,t,dt,xt0,tspgm,fspgm,spgm,SUT,analysisWin,Fs,maxIterationVoted,{'Inconsistency','frogE'},plotIter,filtm,fMaxStftFilter,plotFilt,f);
                        
                        
                    end
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    xt1=xt;%
                    %
                    % figure;
                    % plot(t,real(xt));
                    % hold on ;
                    % plot(t,abs(xt));
                    % yyaxis right;
                    % plot(t,unwrap(angle(xt)))
                    %
                    
                    
                    
                    
                    
                    marginal=sum(abs(Si).^2,1);
                    clear Si S0 iniPhase SiOut SiOut_Raw spgm
                    
                    
                    
                    
                    %% Phase analysis
                    %
                    % %
                    % %
                    % sec1t=118:554;
                    % sec1f=192:257;
                    % % xt=filtSG_tf(xt1,t,f,fMaxStftFilter/df,2,1);
                    % % sec1t=554:1391;
                    % % sec1f=391:465;
                    %
                    % t1=t(sec1t);
                    % xt1=xt(sec1t);
                    % f1=linspace(-Fs/2,Fs/2,numel(t1));
                    % xf1=(nfft(xt1));
                    % fitObj1=fit(2*pi*f1(sec1f)',unwrap(angle(xf1(sec1f)))','poly2');
                    %
                    % figure;plot(2*pi*f1,unwrap(angle(nfft(xt1)))); hold on; plot(2*pi*f1(sec1f),fitObj1(2*pi*f1(sec1f)));
                    % figure;plot(unwrap(angle((xt1))));
                    % beta2=fitObj1.p1*2
                    % beta2/22e-24
                    % %% Phase analysis
                    
                    
                    if phaseAnalysis==1
                        %     xt=xt(4111:9808);t=t(4111:9808);
                        recovPhaseRaw=unwrap(angle(xt));
                        linFit=fit(t',recovPhaseRaw','poly1')
                        recovPhaseflat=recovPhaseRaw-linFit(t)';
                        recovPhase=recovPhaseflat;
                        % recovPhase=real(filtSG_tf(recovPhaseflat,t,f,200,2,1));
                        
                        %     figure;plot(recovPhase)
                        %  fitReg={1:721,860:1590};
                        
                        [~,locsTop]=findpeaks(recovPhase,'MinPeakProminence',50);
                        [~,locsBot]=findpeaks(-recovPhase,'MinPeakProminence',50);
                        
                        rangeTop=round(7*winLen_t/dt); rangeBot=round(3*winLen_t/dt);
                        topInds=((1:rangeTop)-rangeTop/2)+locsTop';
                        botInds=((1:rangeBot)-rangeBot/2)+locsBot';
                        fitReg={topInds,botInds};
                        
                        figure;
                        subplot(3,1,1)
                        plot(t*1e9,recovPhase); hold on
                        xlim([t(1) t(end)]*1e9);
                        
                        for ifit=1:2
                            %  iFit=2;
                            
                            tfits=t(fitReg{ifit}); yfits=recovPhase(fitReg{ifit});
                            beta2=[];
                            for i=1:numel(tfits(:,1))
                                tfit=tfits(i,:); yfit=yfits(i,:);
                                fitObj=fit(tfit',yfit','poly2')
                                
                                beta2(i)=1/(2*fitObj.p1)
                                %  pause(0.3)
                                plot(tfit*1e9,yfit); plot(tfit*1e9,fitObj(tfit));plot( tfit(round(numel(yfit)/2))*1e9,yfit(round(numel(yfit)/2)),'*')
                            end
                            beta2s{ifit}=beta2;
                        end
                        
                        
                        
                        beta2s1=beta2s{1}; beta2s1(isoutlier(beta2s1,'gesd',"MaxNumOutliers",1))=[];
                        beta2s2=beta2s{2}; beta2s2(isoutlier(beta2s2,'gesd',"MaxNumOutliers",1))=[];
                        
                        meanb2(1)=mean(beta2s1)*1e24;
                        stdb2(1)=std(beta2s1)*1e24;
                        meanb2(2)=mean(beta2s2)*1e24;
                        stdb2(2)=std(beta2s2)*1e24;
                        
                        ylabel('Phase (rad)');
                        yyaxis right
                        plot(t*1e9,abs(xt).^2/(max(abs(xt).^2)))
                        xlabel('time (ns)')
                        
                        subplot(3,1,2)
                        plot(t(topInds(:,rangeTop/2))*1e9,beta2s{1},'*')
                        hold on
                        plot(t(topInds(:,rangeTop/2))*1e9,mean(beta2s{1})*ones(1,numel(beta2s{1})))
                        title(['mean beta2: ' num2str(meanb2(1)) '±' num2str(stdb2(1)) ...
                            'ps2 (' num2str(meanb2(1)/22) '±' num2str(stdb2(1)/22) ' km SMF)'])
                        xlim([t(1) t(end)]*1e9);
                        
                        
                        subplot(3,1,3)
                        plot(t(botInds(:,rangeBot/2))*1e9,beta2s{2},'*')
                        hold on
                        plot(t(botInds(:,rangeBot/2))*1e9,mean(beta2s{2})*ones(1,numel(beta2s{2})))
                        title(['mean beta2: ' num2str(meanb2(2)) '±' num2str(stdb2(2)) ...
                            'ps2 (' num2str(meanb2(2)/22) '±' num2str(stdb2(2)/22) ' km SMF)']);
                        xlim([t(1) t(end)]*1e9);
                        
                        
                        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR PARAMETRIC ANALYSIS
                        %     allMeanb2s1(icl, inptWin, ibwIRF, iInter)=meanb2(1);
                        % allMeanb2s2(icl, inptWin, ibwIRF, iInter)=meanb2(2);
                        % allStdb2s1(icl, inptWin, ibwIRF, iInter)=stdb2(1);
                        % allStdb2s2(icl, inptWin, ibwIRF, iInter)=stdb2(2);
                        %     timeNowStr=getTimeString();
                        % itStr=num2str([icl inptWin ibwIRF iInter])
                        % save_figpng([outputFigsDir 'beta2_' itStr '_' timeNowStr]);
                        
                        
                    end
                    
                    
                    
                    if phaseAnalysis==2
                        
                        
                        figure;
                        subplot(3,1,1)
                        plot(t-t(end)/2,abs(xt).^2); ylabel('intensity')
                        hold on
                        plot(tspgm-mean((tspgm)),marginal/max(marginal)*max(abs(xt).^2))
                        yyaxis right
                        plot(t-t(end)/2,unwrap(angle(xt))); ylabel('angle (rad)')
                        xlabel('time (ns)')
                        subplot(3,1,2)
                        
                        movMeanXt=movmean(abs(xt).^2,winLen/2);
                        freqMovMean=abs(nfft(movMeanXt));
                        [~,DCpeak]=max(freqMovMean); freqMovMean([DCpeak-1 DCpeak DCpeak+1])=0;
                        freqMovMeanFilt=smooth(abs(freqMovMean),winLen/2).';
                        [maxFreq]=find(freqMovMeanFilt>max(freqMovMeanFilt)*0.1,1);
                        plot(t,movMeanXt);
                        subplot(3,1,3)
                        plot(f,freqMovMean)
                        hold on
                        plot(f,freqMovMeanFilt)
                        plot(f(maxFreq),freqMovMean(maxFreq),'*')
                        xlim([-2*abs(f(maxFreq)) 2*abs(f(maxFreq))])
                        
                        
                        
                        %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR PARAMETRIC ANALYSIS
                        allMaxFreqs( icl, ifn, inptWin, ibwIRF, iInter)=abs(f(maxFreq));
                        
                        timeNowStr=getTimeString();
                        %                         save_figpng([outputFigsDir 'temporalTrace_' itStr '_' timeNowStr]);
                        
                        f1=gcf;
                        saveas(f1,[outputFigsDir 'temporalTrace_' itStr '_' timeNowStr '.png'],'png')
                        
                        
                        
                    end
                    
                    
                    
                    
                    
                    
                    close all
                end
            end
        end
        
        timeNowStr=getTimeString();
        save([outputFigsDir 'dataMatrix' timeNowStr],'allMaxFreqs','alldiCs','alldiRs')
        %
    end
    
    
    timeNowStr=getTimeString();
    save([outputFigsDir 'dataMatrixAll' timeNowStr],'allMaxFreqs','alldiCs','alldiRs')
    %
end

