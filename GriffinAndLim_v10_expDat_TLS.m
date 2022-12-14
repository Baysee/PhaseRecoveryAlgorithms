% Spectrogram function
addpath( 'C:\Users\Lord Photon\Documents\MATLAB\library_repo\library' )
addpath( '/Users/ben/Documents/MATLAB/library_repo' )
addpath('../phaseRecovery_Data');
outputFigDir='../phaseRecov_Results/'
% save_figpng
%% Load data and set time-frequency vectors


% load('RTOZigZag.mat');
% load('/Users/ben/Documents/MATLAB/timeFrequencyAnalysis/phaseRecovery_Data/RTOzigzag_deconv.mat');
% load('C:\Users\Lord Photon\Documents\MATLAB\time-frequency analysis\PhaseRecoveryAlgorithms_repo\phaseRecovery_Data/RTOZigZag.mat');
% load('RTOzigzag_deconv.mat');
load('RTOzigzag_deconv_20221028.mat');
lowerClip=max(max(spgm))/70;
tIndsExpInterest=2:179;
fSpecGHz=fSpecGHz;%/(56.4e9*TargetResolution)
winLent=200e-12
nptPerWin=32;
phaseAnalysis=1;
nFreqElem=numel(fSpecGHz);%2000;
freqInds=round(linspace(1,numel(fSpecGHz),nFreqElem));
filtm=10;
fMaxStft=80e9;
bwIRF=50e9;
plotFilt=0;

% % % % yr_f=nfft(yReshaped,1); IRF_reconv_filt=superGauss(0,bwIRF,1,fRaw,0);
% 


% % % 




% 
% load('OsoDeconv_narrow_20221028.mat');spgm=circshift(spgm,-82);lowerClip=max(max(spgm))/60;
% load('OSO_deconvWithBackground.mat');spgm=circshift(spgm,-82);lowerClip=max(max(spgm))/60;
% load('OSO_deconvNoBackground.mat');spgm=circshift(spgm,-82);lowerClip=max(max(spgm))/60;
% load('OSOdataCross_filtv2.mat');lowerClip=max(max(spgm))/25;
% load('OSOdataCross.mat');lowerClip=max(max(spgm))/25;
% tIndsExpInterest=40:200;
% winLent=62.5e-12
% bwIRF=500e9;
% nptPerWin=64;
% phaseAnalysis=2;
% fMaxStft=400e9;
% 










dtSpgmRaw=winLent/numel(spgm(:,1));
tRaw=(1:numel(dtSpgmRaw))*dtSpgmRaw;
fRaw=linspace(-1/(2*dtSpgmRaw),1/(2*dtSpgmRaw),numel(spgm));
yReshaped=reshape(spgm,[1,numel(spgm)]);
yReconv=abs(filterSig(1,bwIRF,yReshaped,1,fRaw));
spgm=abs(reshape(yReconv,size(spgm)));



fMaxTLS=fSpecGHz(end)*1e9*2;
spgmIni=spgm;%.^0.8;%(freqInds,:);
% restrict time axis and clip lower values to avoid noise issues.
spgmExpRaw=spgmIni(:,tIndsExpInterest);
spgmExpRaw(spgmExpRaw<lowerClip)=0;%
spgmExpRaw(spgmExpRaw>lowerClip)=spgmExpRaw(spgmExpRaw>lowerClip)-lowerClip;
spgm=spgmExpRaw;

tspgm=tSpecns(tIndsExpInterest)*1e-9;%linspace(tspgm_raw(1),tspgm_raw(end),numel(tspgm_raw)*interpAmount_t);
fspgm=fSpecGHz*1e9;%linspace(fspgm_raw(1),fspgm_raw(end),numel(fspgm_raw)*interpAmount_f);%fspgm_raw;

dt=winLent/numel(spgm(:,1));
%% Time frequency vectors definition

% 
% fSpecExp=fSpecGHz*1e9;  tSpecExp=tSpecns(tIndsExpInterest)*1e-9;
% % Setup time-frequency axes
% %2^nextpow2(numel(fSpecExpRaw));
% %%% Used to have a +1 in the tWind???? % tWind=winLen_t*round(1+(tSpecExp(end)-tSpecExp(1))/winLen_t); % The total length of t dictates the minimum required frequency resolution
% tWind=numel(tSpecExp)*winLen_t;%(tSpecExp(end)-tSpecExp(1));%winLen_t*round((tSpecExp(end)-tSpecExp(1))/winLen_t); % The total length of t dictates the minimum required frequency resolution
% TargetResolution=winLen_t/nptPerWin;%1/(extraFreqNeeded*2*fSpecExpRaw(end));%448e9;%winLen_t/(2^5); % Target temporal resolution of the STFT (i.e., frequency span)
% % % %  Used to force to round by 2 % % % freqPadNeeded=round((1/TargetResolution-2*fSpecExpRaw(end))*tWind/2)*2; %%%%%%%%%% rounding to 2 because I'm lazy now. is this important?
% % freqPadNeeded=round((1/TargetResolution-2*fSpecExpRaw(end))/2*tWind)*2; %%%%%%%%%% rounding to 2 because I'm lazy now. is this important?
% 
% lent=round(tWind/TargetResolution);
% 
% % % 0pad in the time domain
% % nSlicesPad=2^8;
% % tPad=nptPerWin*nSlicesPad;
% % spgmExpRaw=[zeros(numel(fSpecExp),nSlicesPad/2),spgmExpRaw,zeros(numel(fSpecExp),nSlicesPad/2)];
% % lent=tPad+lent; tSpecExp=(1:numel(spgmExpRaw(1,:)))*winLen_t;tWind=tWind+nSlicesPad*winLen_t;
% 
% dt=tWind/lent;tWind=lent*dt;
% t=(1:lent)*dt;
% Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
% fG=f*10^-9;tps=t*10^12;%GHz
% scale=1;
% 
% 
% fspgm_raw=f;
% tspgm_raw=tSpecExp;
% SUT=zeros(1,lent);
% % Zero pad if needed
% if 1/(2*TargetResolution)>fSpecExp(end)
%     
% [a,nZeros]=find(f>-abs(fMaxTLS)/2,1);
% f_tls1=[f(1:nZeros-1),fSpecExp,-flip(f(1:nZeros-1))];
% padded1=[zeros(nZeros-1,numel(spgmExpRaw(1,:)));spgmExpRaw;zeros(nZeros-1,numel(spgmExpRaw(1,:)))];
% spgmExpRawPadded=interp2fun(tspgm_raw,f_tls1,padded1,tspgm_raw,f);
% % spgmExpRawPadded=imresize(padded1,[numel(tspgm_raw),numel(f)]);
% spgmExpRawPadded(1:nZeros,:)=0;spgmExpRawPadded(lent-nZeros:end,:)=0;
% 
% 
% % fSpecExpRawPadded=[-1/(2*TargetResolution),fSpecExp,1/(2*TargetResolution)];%(1:lent)*dfSER-1/(2*TargetResolution);
% % spgmExpRawPadded=[zeros(1,numel(tSpecExp));spgmExpRaw;zeros(1,numel(tSpecExp))]; % If freqPadNeeded/2 doesn't work this will cause an error
% else
% fSpecExpRawPadded=fSpecExp; spgmExpRawPadded=spgmExpRaw;
% end
% 
% 
% 

%% window Setup

winLen=numel(spgm(:,1));
lent=numel(spgm);
winInc=winLen;
interpAmount_t=1 ; % For now, make this a power of 2 (or 1)!!
interpAmount_f=1  ; % For now, make this a power of 2 (or 1)!!


winSec=ones(1,winLen);
% win=superGauss(0,winLen/2,2,windowInds,0);
win=zeros(1,lent); win(round(lent/2)-floor(winLen/2):round(lent/2)+ceil(winLen/2)-1)=winSec;
win=circshift(win,round(lent/2));

% No need to change the ones below
nIncs=lent/winInc; % By making winInc a power of 2, we can make sure to have an integer number of windows.
windowCenters=(1:nIncs)*winInc;



t=(1:lent)*dt;
Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
fG=f*10^-9;tps=t*10^12;%GHz
scale=1;


CL=2*pi*fMaxTLS/winLent;
phi2=1/CL;

phit=repmat( ((1:winLen)/winLen*winLent-winLent/2).^2*CL/2, [1,nIncs]) ;
phiw=phi2/2*(2*pi*f).^2;





% fspgm=fspgm_raw;
% fspgm=fSpecExpRawPadded;
% spgm=griddata(tspgm_rawM,fspgm_rawM,spgmRaw,tspgmM,fspgmM,'natural');
% spgminterp1=interp2fun(tspgm_raw,fspgm,spgmRaw,tspgm,fspgm);
% 
% spgm=spgminterp1;

% % % % spgm=abs(interp2fun(tspgm_raw,fSpecExpRawPadded,spgmExpRawPadded,tspgm,fspgm));
% if interpAmount_t>1
%      spgm=abs(imresize(spgmExpRawPadded,[lent,numel(tspgm)]));
% else
%     spgm=spgmExpRawPadded;
% end
% 
% % Remove the extrapolated noise
% numToDecimate=round((1/(2*dt)-fSpecExp(end))/df);
% spgm(1:numToDecimate,:)=0; spgm(end-numToDecimate:end,:)=0;



%% for smoothing out the spgm
% sigma=1;
% %   lowerClip=7.7;
% % spgm(spgm<lowerClip)=0;%
% % spgm(spgm>lowerClip)=spgm(spgm>lowerClip)-lowerClip;
% spgminterp1=spgminterp1-mean(spgminterp1(:,1));
% spgm= imgaussfilt(spgminterp1,sigma);



figure;
% subplot(3,1,1)
% imagesc(tSpecns,fSpecGHz,spgmIni);
% % subplot(3,1,2)
% % imagesc(spgminterp1)
% subplot(3,1,3)
imagesc(tspgm,fspgm,spgm)





%% Iterative Griffin and Lim algorithm

nIter=3;
SiOut=zeros(numel(spgm(:,1)),numel(spgm(1,:)),nIter);
xtOut=zeros(numel(spgm),nIter);
for iIter=1:nIter

S0=sqrt(spgm).*exp(1j*rand(size(spgm))*2*pi);%.*(-1*(stft<0));%.*exp(1j*rand(size(spgm))*2*pi); % Seed stft

% figure;plot(tspgm,sum(S0,1))

xt0=istft_TLS(S0,winLen,phit,phiw);;%get_istft_fullSigLen(lent,windowCenters,analysisWin,Fs,nIncs,S0);
SUT=xt0;
maxIteration=40;

% Convergence criterion
di=zeros(1,maxIteration);
diC=di;
diR=di;
plotIter=1;

if iIter==1; plotIter=1;end

[xt,Si]=phaseRecovLoopTLS(nIncs,windowCenters,lent,winLen,t,f,dt,xt0,tspgm,fspgm,spgm,SUT,Fs,maxIteration,phit,phiw);




SiOut(:,:,iIter)=Si;
xtOut(:,iIter)=xt;
iIter
end



SiOut_Raw=SiOut;









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
% xt0=get_istft_fullSigLen(lent,windowCenters,analysisWin,Fs,nIncs,SiPairOpt);
xt0=istft_TLS(SiPairOpt,winLen,phit,phiw);;%get_istft_fullSigLen(lent,windowCenters,analysisWin,Fs,nIncs,S0);

plotIter=1;
% xt0=mean(xtOut,2).';
maxIteration=40;
% xt0=get_istft_fullSigLen(lent,windowCentersInterp,analysisWin,Fs,nIncsInterp,sqrt(spgm).*exp(1j*angle(SiMean)));
[xt,Si]=phaseRecovLoopTLS(nIncs,windowCenters,lent,winLen,t,f,dt,xt0,tspgm,fspgm,spgm,SUT,Fs,maxIteration,phit,phiw);















xt1=xt;%
%
% figure;
% plot(t,real(xt));
% hold on ;
% plot(t,abs(xt));
% yyaxis right;
% plot(t,unwrap(angle(xt)))
%










%% Phase analysis 

%
%
sec1t=118:554;
sec1f=192:257;
xt=filtSG_tf(xt1,t,f,fMaxTLS/df,2,1);
% sec1t=554:1391;
% sec1f=391:465;

t1=t(sec1t);
xt1=xt(sec1t);
f1=linspace(-Fs/2,Fs/2,numel(t1));
xf1=(nfft(xt1));
fitObj1=fit(2*pi*f1(sec1f)',unwrap(angle(xf1(sec1f)))','poly2');

figure;plot(2*pi*f1,unwrap(angle(nfft(xt1)))); hold on; plot(2*pi*f1(sec1f),fitObj1(2*pi*f1(sec1f)));
figure;plot(unwrap(angle((xt1))));
beta2=fitObj1.p1*2
beta2/22e-24
%% Phase analysis


if phaseAnalysis==1
%     xt=xt(4111:9808);t=t(4111:9808); 
recovPhaseRaw=unwrap(angle(xt));
linFit=fit(t',recovPhaseRaw','poly1')
recovPhaseflat=recovPhaseRaw-linFit(t)';
recovPhase=recovPhaseflat;
% recovPhase=real(filtSG_tf(recovPhaseflat,t,f,200,2,1));

 figure;plot(recovPhase)
%  fitReg={1:721,860:1590};

 [~,locsTop]=findpeaks(recovPhase,'MinPeakProminence',50);
 [~,locsBot]=findpeaks(-recovPhase,'MinPeakProminence',50);
 
 rangeTop=round(7*winLent/dt); rangeBot=round(3*winLent/dt);
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
 ylabel('Phase (rad)');
yyaxis right
plot(t*1e9,abs(xt).^2/(max(abs(xt).^2)))
 xlabel('time (ns)')

  subplot(3,1,2)
plot(t(topInds(:,rangeTop/2))*1e9,beta2s{1},'*')
hold on
plot(t(topInds(:,rangeTop/2))*1e9,mean(beta2s{1})*ones(1,numel(beta2s{1})))
title(['mean beta2: ' num2str(mean(beta2s{1})*1e24) 'ps2 (' num2str(mean(beta2s{1})*1e24/22) ' km SMF)'])
 xlim([t(1) t(end)]*1e9);

   subplot(3,1,3)
plot(t(botInds(:,rangeBot/2))*1e9,beta2s{2},'*')
hold on
plot(t(botInds(:,rangeBot/2))*1e9,mean(beta2s{2})*ones(1,numel(beta2s{2})))
title(['mean beta2: ' num2str(mean(beta2s{2})*1e24) 'ps2 (' num2str(mean(beta2s{2})*1e24/22) ' km SMF)'])
 xlim([t(1) t(end)]*1e9);

end



if phaseAnalysis==2
    
    
    figure;
    plot(t,abs(xt).^2); ylabel('intensity')
    yyaxis right
    plot(t,unwrap(angle(xt))); ylabel('angle (rad)')
    xlabel('time (ns)')
    
end
