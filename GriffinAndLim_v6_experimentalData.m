% Spectrogram function
addpath( 'C:\Users\Lord Photon\Documents\MATLAB\library_repo\library' )
addpath( '/Users/ben/Documents/MATLAB/library_repo' )
%% Load data and set time-frequency vectors

% % load('/Users/ben/Documents/MATLAB/timeFrequencyAnalysis/phaseRecovery_Data/OSOdataCross_filtv2.mat');
% load('/Users/ben/Documents/MATLAB/timeFrequencyAnalysis/phaseRecovery_Data/RTOZigZag.mat');
% load('C:\Users\Lord Photon\Documents\MATLAB\time-frequency analysis\PhaseRecoveryAlgorithms_repo\phaseRecovery_Data/RTOZigZag.mat');
load('C:\Users\Lord Photon\Documents\MATLAB\time-frequency analysis\PhaseRecoveryAlgorithms_repo\phaseRecovery_Data/RTOzigzag_deconv.mat');
% tIndsExpInterest=68:129;
tIndsExpInterest=1:186;
fSpecGHz=fSpecGHz;%/56.4*55;
winLen_t=200e-12
lowerClip=0.00000;
extraFreqNeeded=1;
phaseAnalysis=1;

% 
% load('/Users/ben/Documents/MATLAB/timeFrequencyAnalysis/phaseRecovery_Data/OSOdataCross.mat');
% tIndsExpInterest=70:160;
% winLen_t=62.5e-12
% lowerClip=8;
% extraFreqNeeded=1;
% phaseAnalysis=2;




spgmIni=spgm.^2;
% restrict time axis and clip lower values to avoid noise issues.

spgmExpRaw=spgmIni(:,tIndsExpInterest);
spgmExpRaw(spgmExpRaw<lowerClip)=0;%
spgmExpRaw(spgmExpRaw>lowerClip)=spgmExpRaw(spgmExpRaw>lowerClip)-lowerClip;


fSpecExp=fSpecGHz*1e9;  tSpecExp=tSpecns(tIndsExpInterest)*1e-9;
% Setup time-frequency axes
%2^nextpow2(numel(fSpecExpRaw));
%%% Used to have a +1 in the tWind???? % tWind=winLen_t*round(1+(tSpecExp(end)-tSpecExp(1))/winLen_t); % The total length of t dictates the minimum required frequency resolution
tWind=numel(tSpecExp)*winLen_t;%(tSpecExp(end)-tSpecExp(1));%winLen_t*round((tSpecExp(end)-tSpecExp(1))/winLen_t); % The total length of t dictates the minimum required frequency resolution
TargetResolution=10e-12;%1/(extraFreqNeeded*2*fSpecExpRaw(end));%448e9;%winLen_t/(2^5); % Target temporal resolution of the STFT (i.e., frequency span)
% % %  Used to force to round by 2 % % % freqPadNeeded=round((1/TargetResolution-2*fSpecExpRaw(end))*tWind/2)*2; %%%%%%%%%% rounding to 2 because I'm lazy now. is this important?
% freqPadNeeded=round((1/TargetResolution-2*fSpecExpRaw(end))/2*tWind)*2; %%%%%%%%%% rounding to 2 because I'm lazy now. is this important?

lent=round(tWind/TargetResolution/2)*2;
dt=tWind/lent;
dfSER=1/tWind;%;fSpecExpRaw(2)-fSpecExpRaw(1);
fSpecExpRawPadded=[-1/(2*TargetResolution),fSpecExp,1/(2*TargetResolution)];%(1:lent)*dfSER-1/(2*TargetResolution);
% fSpexExpRawPadded=([1:freqPadNeeded+numel(fSpecExpRaw)]-round((freqPadNeeded+numel(fSpecExpRaw))/2))*dfSER;

% spgmExpRawPadded=[zeros(freqPadNeeded/2,numel(tSpecExp));spgmExpRaw;zeros(freqPadNeeded/2,numel(tSpecExp))]; % If freqPadNeeded/2 doesn't work this will cause an error
spgmExpRawPadded=[zeros(1,numel(tSpecExp));spgmExpRaw;zeros(1,numel(tSpecExp))]; % If freqPadNeeded/2 doesn't work this will cause an error


% lent=round(tWind/TargetResolution);
t=(1:lent)*dt;
Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
fG=f*10^-9;tps=t*10^12;%GHz
scale=1;


fspgm_raw=f;
tspgm_raw=tSpecExp;
nWindsNoOverlap=numel(tSpecExp);
% spgmRaw=interp2fun(tSpecExp,fSpecExpRawPadded,spgmExpRawPadded,tspgm_raw,fspgm_raw);
spgmRaw=interp2fun(tSpecExp,fSpecExpRawPadded,spgmExpRawPadded,tspgm_raw,fspgm_raw);

numToDecimate=round((1/(2*dt)-fSpecExp(end))/df);
spgmRaw(1:numToDecimate,:)=0; spgmRaw(end-numToDecimate:end,:)=0;

%% SUT generation
% Used to compared with a simulation SUT... doesn't seem very relevant.

fmax=190e9;
SUTf=superGauss(0,fmax,4,f,0).*(exp(1j*(105*22e-24)*(2*pi*f).^2/2))+superGauss(0,fmax,10,f,0).*(exp(-1j*(105*22e-24)*(2*pi*f).^2/2));
% SUTf=superGauss(0,sutBW,10,f,0).*(exp(1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2))+...
%     superGauss(0,sutBW,10,f,0).*(exp(-1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2));
SUT=nifft(SUTf,Fs);


% spgmRaw=abs(get_stft_fullSigLen(nIncs,winInds,windowCenters,lent,win,dt,winLen,SUT));
% spgmRaw=abs(get_stft_fullSigLen(nIncs,windowCenters,lent,win,dt,SUT));




%% window Setup

% Adjust these parameters as needed
winLen=round(winLen_t/dt);
if (winLen_t/dt-round(winLen_t/dt))>0.01
    warning('big rounding off winLen_t?')
end
winInc=winLen;%winLen-1;%/(2^2);
interpAmount_t=1  ; % For now, make this a power of 2 (or 1)!!
interpAmount_f=1; % For now, make this a power of 2 (or 1)!!

% win=hann(winLen+2).^3;win=win(2:end-1)';%ones(1,winLen);

winSec=ones(1,winLen);

windowInds=(1:lent)-lent/2;
% win=superGauss(0,winLen/2,2,windowInds,0);

win=zeros(1,lent); win(round(lent/2)-winLen/2:round(lent/2)+winLen/2-1)=winSec;
win=circshift(win,lent/2);

% No need to change the ones below
nIncs=lent/winInc; % By making winInc a power of 2, we can make sure to have an integer number of windows.
windowCenters=(1:nIncs)*winInc;
winInds=(1:winLen)-winLen/2;




%% Setup interpolation
nIncsInterp=nIncs*interpAmount_t;
windowCentersInterp=(1:nIncsInterp)*winInc/interpAmount_t;

tspgm=linspace(tspgm_raw(1),tspgm_raw(end),numel(tspgm_raw)*interpAmount_t);
fspgm=linspace(fspgm_raw(1),fspgm_raw(end),numel(fspgm_raw)*interpAmount_f);%fspgm_raw;

% spgm=griddata(tspgm_rawM,fspgm_rawM,spgmRaw,tspgmM,fspgmM,'natural');
spgminterp1=interp2fun(tspgm_raw,fspgm,spgmRaw,tspgm,fspgm);

spgm=spgminterp1;
% 
% sigma=1;
% %   lowerClip=7.7;
% % spgm(spgm<lowerClip)=0;%
% % spgm(spgm>lowerClip)=spgm(spgm>lowerClip)-lowerClip;
% spgminterp1=spgminterp1-mean(spgminterp1(:,1));
% spgm= imgaussfilt(spgminterp1,sigma);



figure;
subplot(3,1,1)
imagesc(spgmIni);
subplot(3,1,2)
imagesc(spgminterp1)
subplot(3,1,3)
imagesc(tspgm,fspgm,spgm)



winLenInterp=numel(fspgm)*interpAmount_f;
winInterp=interp1(linspace(0,1,lent),win,linspace(0,1,winLenInterp));
winIndsInterp=(1:numel(fspgm))-round(numel(fspgm)/2);
overlapAmount=interpAmount_t*(sum(win))/winInc;%numel(winIndsInterp)/(windowCentersInterp(2)-windowCentersInterp(1)); % This is the "overlapamount" AFTER interpolation
analysisWin=winInterp/overlapAmount; % Analysis window for the inverse spgm


%% Iterative Griffin and Lim algorithm

S0=sqrt(spgm);%.*exp(1j*rand(size(spgm))*2*pi);%.*(-1*(stft<0));%.*exp(1j*rand(size(spgm))*2*pi); % Seed stft

xt0=get_istft_fullSigLen(lent,windowCentersInterp,analysisWin,Fs,nIncsInterp,S0);

maxIteration=30;

% Convergence criterion
di=zeros(1,maxIteration);
diC=di;
diR=di;
[xt,Si]=phaseRecovLoop(nIncsInterp,windowCentersInterp,lent,winInterp,winLen,t,dt,xt0,tspgm,fspgm,spgm,SUT,analysisWin,Fs,maxIteration);
%
%
% figure;
% plot(t,real(xt));
% hold on ;
% plot(t,abs(xt));
% yyaxis right;
% plot(t,unwrap(angle(xt)))
%
%
%



%% Phase analysis


if phaseAnalysis==1
    
recovPhaseRaw=unwrap(angle(xt));
linFit=fit(t',recovPhaseRaw','poly1')
recovPhaseflat=recovPhaseRaw-linFit(t)';
recovPhase=real(filtSG_tf(recovPhaseflat,t,f,200,2,1));

 figure;plot(recovPhase)
%  fitReg={1:721,860:1590};

 [~,locsTop]=findpeaks(recovPhase,'MinPeakProminence',50);
 [~,locsBot]=findpeaks(-recovPhase,'MinPeakProminence',50);
 
 rangeTop=200; rangeBot=100;
 topInds=((1:rangeTop)-rangeTop/2)+locsTop';
 botInds=((1:rangeBot)-rangeBot/2)+locsBot';
 fitReg={topInds,botInds};
 
 
 figure;
 subplot(3,1,1)
 plot(t,recovPhase); hold on
 for ifit=1:2
%  iFit=2;
 
 tfits=t(fitReg{ifit}); yfits=recovPhase(fitReg{ifit});
beta2=[];
 for i=1:numel(tfits(:,1))
     tfit=tfits(i,:); yfit=yfits(i,:);
 fitObj=fit(tfit',yfit','poly2')
 
beta2(i)=2/fitObj.p1
plot(tfit,yfit); plot(tfit,fitObj(tfit));plot( tfit(round(end/2)),yfit(round(end/2)),'*')
 end
 beta2s{ifit}=beta2;
 end
yyaxis right
plot(t,abs(xt).^2)
 
 
  subplot(3,1,2)
plot(topInds(:,rangeTop/2),beta2s{1},'*')
hold on
plot(topInds(:,rangeTop/2),mean(beta2s{1})*ones(1,numel(beta2s{1})))
title(['mean beta2: ' num2str(mean(beta2s{1})*1e24) 'ps2 (' num2str(mean(beta2s{1})*1e24/22) ' km SMF)'])

   subplot(3,1,3)
plot(botInds(:,rangeBot/2),beta2s{2},'*')
hold on
plot(botInds(:,rangeBot/2),mean(beta2s{2})*ones(1,numel(beta2s{2})))
title(['mean beta2: ' num2str(mean(beta2s{2})*1e24) 'ps2 (' num2str(mean(beta2s{2})*1e24/22) ' km SMF)'])

end



if phaseAnalysis==2
    
end
