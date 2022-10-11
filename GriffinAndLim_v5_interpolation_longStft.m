% Spectrogram function
 addpath( 'C:\Users\Lord Photon\Documents\MATLAB\library_repo\library' )
addpath( '/Users/ben/Documents/MATLAB/library_repo' )
%% Time frequency vectors definition

lent=2^13;                      % Signal length
tWind=4e-9;                   % Time window span


t=linspace(0,tWind,lent);
dt=t(2)-t(1);Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
fG=f*10^-9;tps=t*10^12;%GHz
scale=1;



%% SUT generation

fmax=31e9/2;%Fs/10;
% SUTf=superGauss(0,fmax,10,f,0).*(exp(1j*(tWind/4/(fmax*2*pi))*(2*pi*f).^2/2));
SUTf=superGauss(0,fmax,10,f,0).*(exp(1j*(240*22e-24/2)*(2*pi*f).^2/2));%+superGauss(0,fmax,10,f,0).*(exp(-1j*(120*22e-24/2)*(2*pi*f).^2/2));
% SUTf=superGauss(0,sutBW,10,f,0).*(exp(1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2))+...
%     superGauss(0,sutBW,10,f,0).*(exp(-1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2));
SUT=nifft(SUTf,Fs);

% SUT=ones(size(SUT));



%% window Setup

% Adjust these parameters as needed
winLen=2^7;
winLent=winLen*dt
winInc=winLen/4;%winLen-1;%/(2^2);
interpAmount_t=1; % For now, make this a power of 2 (or 1)!!
interpAmount_f=1; % For now, make this a power of 2 (or 1)!!


winSec=ones(1,winLen);
windowInds=(1:lent)-lent/2;
% win=superGauss(0,winLen/2,100,windowInds,0);
win=zeros(1,lent); win(round(lent/2)-winLen/2:round(lent/2)+winLen/2-1)=winSec;

win=circshift(win,lent/2);


% No need to change the ones below
nIncs=round(lent/winInc); % By making winInc a power of 2, we can make sure to have an integer number of windows.
windowCenters=(1:nIncs)*winInc;






%% Spectrogram Algorithm

% Get spgm from windowIncrease Above
stft=get_stft_fullSigLen(nIncs,windowCenters,lent,win,dt,SUT);
spgmRaw=abs(stft).^2;
fspgm_raw=f;
tspgm_raw=linspace(t(1),t(end),numel(stft(1,:)));


%%% PLaying around with COLA
% sutRecon=get_istft_fullSigLen(lent,windowCenters,win/(sum(win)/winInc),Fs,nIncs,stft);
% figure;
% for i=1:winInc
% plot(circshift(win,windowCenters(i)),'--')
% hold on
% end
% figure;imagesc(spgmRaw);
% figure;plot(sutRecon); hold on; plot(circshift(win,lent/2))
% Setup interpolation

%% Setup interpolation

nIncsInterp=nIncs*interpAmount_t;
windowCentersInterp=(1:nIncsInterp)*winInc/interpAmount_t;

tspgm=linspace(t(1),t(end),numel(stft(1,:))*interpAmount_t);
fspgm=linspace(f(1),f(end),numel(fspgm_raw)*interpAmount_f);%fspgm_raw;

 spgm=interp2fun(tspgm_raw,fspgm,spgmRaw,tspgm,fspgm);


winLenInterp=numel(fspgm)*interpAmount_f;
winInterp=interp1(linspace(0,1,lent),win,linspace(0,1,winLenInterp));
winIndsInterp=(1:numel(fspgm))-round(numel(fspgm)/2);
overlapAmount=interpAmount_t*(sum(win))/winInc;%numel(winIndsInterp)/(windowCentersInterp(2)-windowCentersInterp(1)); % This is the "overlapamount" AFTER interpolation
analysisWin=winInterp/overlapAmount; % Analysis window for the inverse spgm


%% Iterative Griffin and Lim algorithm

S0=sqrt(spgm);%.*exp(1j*rand(size(spgm))*2*pi);%.*(-1*(stft<0));%.*exp(1j*rand(size(spgm))*2*pi); % Seed stft

xt0=get_istft_fullSigLen(lent,windowCentersInterp,analysisWin,Fs,nIncsInterp,S0);

maxIteration=150;

% Convergence criterion
di=zeros(1,maxIteration);
diC=di;
diR=di;
[xt,Si]=phaseRecovLoop(nIncsInterp,windowCentersInterp,lent,winInterp,winLen,t,dt,xt0,tspgm,fspgm,spgm,SUT,analysisWin,Fs,maxIteration);
%


