% Spectrogram function
 addpath( 'C:\Users\Lord Photon\Documents\MATLAB\library_repo\library' )
addpath( '/Users/ben/Documents/MATLAB/library_repo' )
%% Time frequency vectors definition

lent=2^13;                      % Signal length
tWind=6e-9;                   % Time window span


t=linspace(0,tWind,lent);
dt=t(2)-t(1);Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
fG=f*10^-9;tps=t*10^12;%GHz
scale=1;



%% stft parameters

% Adjust these parameters as needed
winLen=2^8;
winLent=winLen*dt
winInc=winLen/2;%winLen-1;%/(2^2);
interpAmount_t=2; % For now, make this a power of 2 (or 1)!!
interpAmount_f=1; % For now, make this a power of 2 (or 1)!!

% win=hann(winLen+2).^3;win=win(2:end-1)';%ones(1,winLen);

winSec=ones(1,winLen);

windowInds=(1:winLen)-winLen/2;
windowInds=(1:lent)-lent/2;
% win=superGauss(0,winLen/2,100,windowInds,0);
% win=zeros(1,lent); win(round(lent/2)-winLen/2:round(lent/2)+winLen/2-1)=winSec;

% windowInds=(1:lent)-lent/2;
win=superGauss(0,winLen/2,200,windowInds,0);
win=circshift(win,lent/2);

% No need to change the ones below
nIncs=round(lent/winInc); % By making winInc a power of 2, we can make sure to have an integer number of windows.
windowCenters=(1:nIncs)*winInc;
winInds=(1:winLen)-winLen/2;


%% SUT generation

fmax=411e9/2;%Fs/10;
% SUTf=superGauss(0,fmax,10,f,0).*(exp(1j*(tWind/4/(fmax*2*pi))*(2*pi*f).^2/2));
SUTf=superGauss(0,fmax,10,f,0).*(exp(1j*(120*22e-24/2)*(2*pi*f).^2/2))+superGauss(0,fmax,10,f,0).*(exp(-1j*(120*22e-24/2)*(2*pi*f).^2/2));
% SUTf=superGauss(0,sutBW,10,f,0).*(exp(1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2))+...
%     superGauss(0,sutBW,10,f,0).*(exp(-1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2));
SUT=nifft(SUTf,Fs);

SUT=ones(size(SUT));


%% Spectrogram Algorithm

% Get spgm from windowIncrease Above
stft=get_stft_fullSigLen(nIncs,windowCenters,lent,win,dt,SUT);
spgmRaw=abs(stft).^2;
fspgm_raw=f;%((1:winLen)-winLen/2)/winLen*Fs;
tspgm_raw=linspace(t(1),t(end),numel(stft(1,:)));


%%% PLaying around with COLA
sutRecon=get_istft_fullSigLen(lent,windowCenters,win/(sum(win)/winInc),Fs,nIncs,stft);
figure;
% for i=1:winInc
% plot(circshift(win,windowCenters(i)),'--')
% hold on
% end
% figure;imagesc(spgmRaw);
figure;plot(sutRecon); hold on; plot(circshift(win,lent/2))
% Setup interpolation
nIncsInterp=nIncs*interpAmount_t;
windowCentersInterp=(1:nIncsInterp)*winInc/interpAmount_t;

tspgm=linspace(t(1),t(end),numel(stft(1,:))*interpAmount_t);
fspgm=linspace(f(1),f(end),numel(fspgm_raw)*interpAmount_f);%fspgm_raw;

% [tspgm_rawM,fspgm_rawM]=meshgrid(tspgm_raw,fspgm);
% [tspgmM,fspgmM]=meshgrid(tspgm,fspgm);
% 
% 
% spgm=interp2(tspgm_rawM,fspgm_rawM,spgmRaw,tspgmM,fspgmM,'spline');

 spgm=interp2fun(tspgm_raw,fspgm,spgmRaw,tspgm,fspgm);
% spgm=griddata(tspgm_rawM,fspgm_rawM,spgmRaw,tspgmM,fspgmM,'natural');


% spgm1=griddedInterpolant(tspgm_rawM,fspgm_rawM,spgmRaw,tspgmM,fspgmM,'linear');
% spgmInterpolant=griddedInterpolant({tspgm_raw,fspgm_raw},spgmRaw','nearest');

% spgmInterp=spgmInterpolant({tspgm,fspgm})';
% spgm=griddata(tspgm_rawM,fspgm_rawM,ab,tspgmM,fspgmM);
% 
% spgm=spgmRaw;%imgaussfilt(spgmInterp,interpAmount_t/2);
% figure;subplot(2,1,1)
% imagesc(spgmRaw)
% subplot(2,1,2)
% imagesc(spgm)

%% Iterative Griffin and Lim algorithm


% ispgm=get_istft(lent,winInds,windowCenters,analysisWin,Fs,nIncs,S);
%
% figure;plot(real(ispgm)); hold on; %plot(imag(ispgm)); %plot(abs(ispgm));
% plot(SUT)
% legend('Real ispgm','real SUT')
winLenInterp=numel(fspgm)*interpAmount_f;
winInterp=interp1(linspace(0,1,lent),win,linspace(0,1,winLenInterp));
winIndsInterp=(1:numel(fspgm))-round(numel(fspgm)/2);
overlapAmount=interpAmount_t*(sum(win))/winInc;%numel(winIndsInterp)/(windowCentersInterp(2)-windowCentersInterp(1)); % This is the "overlapamount" AFTER interpolation
analysisWin=winInterp/overlapAmount; % Analysis window for the inverse spgm

% 
% 
% 
% checkCola_2D=zeros(1,lent);%winLen,overlapAmount*2);
% % figure; hold on
% for i=1:nIncs
% %     checkCola_2D(:,i)=circshift(win,(i-1)*winInc);
%  checkCola_2D=checkCola_2D+circshift(analysisWin,(i-1)*winInc).*circshift(winInterp,(i-1)*winInc);
% % plot(circshift(analysisWin,(i-1)*winInc)); hold on; plot(circshift(circshift(winInterp,lent/2),(i-1)*winInc)); drawnow()
% end
% figure;plot(checkCola_2D)
% 


S0=sqrt(spgm).*exp(1j*rand(size(spgm))*2*pi);%.*(-1*(stft<0));%.*exp(1j*rand(size(spgm))*2*pi); % Seed stft
S0Mag=abs(sqrt(spgm)); % S0 has the correct amplitude, but not the correct phase
% % % % S0=sggm;
xt=get_istft_fullSigLen(lent,windowCentersInterp,analysisWin,Fs,nIncsInterp,S0);
xt0=xt; % This is the first initial guess

maxIteration=100;
i=1;

% Convergence criterion
di=zeros(1,maxIteration);
diC=di;
diR=di;

h1=figure;
h1.Position=[55 117 990 861];%[-1528 112 821 865];
xlimsZoom=t(round(lent/2))+winLen*dt*[-2 2];

smoothingAmount=40;
while i<maxIteration+1
    
    Si=get_stft_fullSigLen(nIncsInterp,windowCentersInterp,lent,winInterp,dt,xt);         % Get spgm of present signal guess xt
    Sip1=S0Mag.*Si./abs(Si);%.*exp(1j*angle(Si));%.*Si./abs(Si);                           % Enforce magnitude along with calculated phase from Si
    Sip1(isnan(Sip1))=0;
    xt=get_istft_fullSigLen(lent,windowCentersInterp,analysisWin,Fs,nIncsInterp,Sip1);
    xt=(smooth(abs(xt),40).').*exp(1j*angle(xt));
    % di(i)=sqrt(sum(sum(abs(abs(Si)-abs(Sip1)).^2))
    % di(i)=norm(abs(abs(Si)-abs(Sip1)).^2)/norm(abs(Si).^2)
%     di(i)=(norm(abs(Si)-abs(Sip1))).^2;%)/norm(abs(Si).^2)
    diC(i)=sqrt(    sum(sum( abs(  abs(S0) - abs(Si)  ).^2))  /    sum(sum(  abs(S0).^2 )) );%)/norm(abs(Si).^2)
%     diC2(i)=sqrt(    sum(sum( abs(  sqrt(abs(S0).^2) - sqrt(abs(Si).^2)  ).^2 ))  /    sum(sum(  abs(S0).^2 )) );%)/norm(abs(Si).^2)
%     diC(i)=sqrt(    norm( abs(  sqrt(abs(S0).^2) - sqrt(abs(Si).^2)  ).^2, 'fro')  /    norm(  abs(S0).^2 , 'fro') );%)/norm(abs(Si).^2)
%     diR(i)=sqrt(sum(abs(SUT-xt).^2)/sum(abs(SUT).^2));%)/norm(abs(Si).^2)
%     xc=max(abs(xcorr(SUT,xt,50)))/(sqrt(sum(abs(SUT).^2)*sum(abs(xt).^2)));
    diR(i)=0;%(1-xc);%)/norm(abs(Si).^2)
    
    % plot(t,real(xt)); drawnow();

    
    if mod(i,5)==1
    updatePlot(h1,xlimsZoom,t,tspgm,fspgm,xt0,SUT,xt,spgm,Si,diC,diR)
    end
    
        i=i+1;
end


    updatePlot(h1,xlimsZoom,t,tspgm,fspgm,xt0,SUT,xt,spgm,Si,diC,diR)











function updatePlot(h1,xlimsZoom,t,tspgm,fspgm,xt0,SUT,xt,spgm,Si,diC,diR)
figure(h1);
FS=16;
subplot(3,2,1)
imagesc(tspgm,fspgm,spgm)
title('SUT spgm')
subplot(3,2,2)
imagesc(tspgm,fspgm,abs(Si).^2)
title('current Iteration spgm')
subplot(3,2,3)
plotIniFin(t,xt0,xt,SUT,FS)
subplot(3,2,4)
plotIniFin(t,xt0,xt,SUT,FS)
ylims=ylim(); xlim(xlimsZoom); ylim(ylims)

subplot(3,2,5:6)
yyaxis left
plot(20*log10(diC))
ylabel('Normalized inconsistency (dB)')
yyaxis right
% plot(20*log10(diC)); hold on
plot(10*log10(diR))
ylabel('10*log10(1-[corr. coeff])')
legend('Convergence Curve (GL)','spectral convergence');%,'reconstruction error')
xlabel('Iteration'); 
set(gca,'FontSize',FS)
drawnow();
%
end




% 
% function ispgm=get_istft(lent,winInds,windowCenters,analysisWin,Fs,nIncs,S)
% % Reconstruct temporal waveform from spgm by the overlap and add technique
% 
% ispgm=zeros(1,lent); % Initialize variable
% 
% % figure;
% for i=1:nIncs
%     
%     ift=nifft(S(:,i),Fs);
% %     sutInds=mod(winInds+windowCenters(i),lent)+1;
%     ispgm=ispgm+circshift(analysisWin,windowCenters(i)).*ift.';%,;
%     
%     
%     
% %     yyaxis left; hold off
% % plot(abs(ispgm)); hold off
% %     yyaxis right; 
% %     plot(abs(ift)); hold on; 
% %     plot(analysisWin/max(analysisWin)*max(abs(abs(nifft(S(:,i),Fs)))))
% %     hold off
% %     title(num2str(i))
% %     drawnow()
% end
% 
% % ispgm=circshift(ispgm,lent/2);
% % ispgm=ispgm;
% end
% 
% 
% 
% 
% 
% 
% function Aq=interp2fun(x,y,A,xq,yq)
% 
% [xm,ym]=meshgrid(x,y);
% [xqm,yqm]=meshgrid(xq,yq);
% Aq=interp2(xm,ym,A,xqm,yqm,'spline');
% 
% 
% 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 








%% PLOTTING 


%t(windowCenters);

% 
% %% Spectrogram Plot
% h0=figure;
% h0.Position=[-1402 147 669 830];
% xlimsZoom=t(round(end/2))+winLen*dt*[-0.5 0.5];
% FS=16;
% 
% subplot(2,2,1)
% plot(t,abs(SUT).^2);
% ylabel('Intensity');
% yyaxis right
% plot(t,angle(SUT));
% xlabel('Time (a.u.)'); ylabel('phase (rad)')
% set(gca,'FontSize',FS)
% 
% subplot(2,2,2)
% plot(t,abs(SUT).^2);
% ylabel('Intensity');
% yyaxis right
% plot(t,angle(SUT));
% xlabel('Time (a.u.)'); ylabel('phase (rad)')
% xlim(xlimsZoom);
% set(gca,'FontSize',FS)
% 
% subplot(2,2,3:4)
% imagesc(tspgm,fspgm,spgm);
% ylim([-fmax fmax])
% xlabel('Time'); ylabel('Frequency');
% set(gca,'FontSize',FS)






