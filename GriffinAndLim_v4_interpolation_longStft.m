% Spectrogram function
 addpath( 'C:\Users\Lord Photon\Documents\MATLAB\library_repo\library' )
addpath( '/Users/ben/Documents/MATLAB/library_repo' )
%% Time frequency vectors definition

lent=2^14;                      % Signal length
tWind=800e-9;                   % Time window span


t=linspace(0,tWind,lent);
dt=t(2)-t(1);Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
fG=f*10^-9;tps=t*10^12;%GHz
scale=1;



%% spgm parameters

% Adjust these parameters as needed
winLen=2^9;
winInc=winLen/2^4;%winLen-1;%/(2^2);
interpAmount_t=1; % For now, make this a power of 2 (or 1)!!
interpAmount_f=1; % For now, make this a power of 2 (or 1)!!

% win=hann(winLen+2).^3;win=win(2:end-1)';%ones(1,winLen);

winSec=ones(1,winLen);

% windowInds=(1:winLen)-winLen/2;
% win=superGauss(0,winLen/2,100,windowInds,0);

win=zeros(1,lent); win(round(lent/2)-winLen/2:round(lent/2)+winLen/2-1)=winSec;
win=circshift(win,lent/2);

% No need to change the ones below
nIncs=lent/winInc; % By making winInc a power of 2, we can make sure to have an integer number of windows.
windowCenters=(1:nIncs)*winInc;
winInds=(1:winLen)-winLen/2;


%% SUT generation

% SUT=exp(1j*pi*sin(2*pi/(40*winLen*dt)*t));



% % % % % % Linearly chirped signal
% fmax=Fs;
% fini=0; ffin=fmax/20;
% c=(ffin-fini)/tWind;
% SUT=exp(1j*pi*sin(2*pi*(c/2*t.^2+fini*t))).*superGauss(0,tWind/3,10,t,tWind/2);
% SUT=sin(2*pi*(c/2*t.^2+fini*t))+1;%.*superGauss(0,tWind/3,10,t,tWind/2);

fmax=Fs/30;
SUTf=superGauss(0,fmax,10,f,0).*(exp(1j*(tWind/4/(fmax*2*pi))*(2*pi*f).^2/2));
% SUTf=superGauss(0,sutBW,10,f,0).*(exp(1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2))+...
%     superGauss(0,sutBW,10,f,0).*(exp(-1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2));
SUT=nifft(SUTf,Fs);
% SUT=superGauss(0,tWind/20,1,t,2*tWind/3)+superGauss(0,tWind/10,2,t,tWind/3);

% SUT=superGauss(0,tWind/3,10,t,tWind/2);




%% Spectrogram Algorithm

% Get stft from windowIncrease Above
spgm=get_spgm(nIncs,winInds,windowCenters,lent,win,dt,winLen,SUT);
stftRaw=abs(spgm).^2;
fstft_raw=f;%((1:winLen)-winLen/2)/winLen*Fs;
tstft_raw=linspace(t(1),t(end),numel(spgm(1,:)));

figure;imagesc(stftRaw);

% Setup interpolation
nIncsInterp=nIncs*interpAmount_t;
windowCentersInterp=(1:nIncsInterp)*winInc/interpAmount_t;

tstft=linspace(t(1),t(end),numel(spgm(1,:))*interpAmount_t);
fstft=linspace(f(1),f(end),numel(fstft_raw)*interpAmount_f);%fstft_raw;

% [tstft_rawM,fstft_rawM]=meshgrid(tstft_raw,fstft);
% [tstftM,fstftM]=meshgrid(tstft,fstft);


% stft=interp2(tstft_rawM,fstft_rawM,stftRaw,tstftM,fstftM,'spline');
% stft=griddata(tstft_rawM,fstft_rawM,stftRaw,tstftM,fstftM,'natural');


% stft1=griddedInterpolant(tstft_rawM,fstft_rawM,stftRaw,tstftM,fstftM,'linear');
% stftInterpolant=griddedInterpolant({tstft_raw,fstft_raw},stftRaw','nearest');

% stftInterp=stftInterpolant({tstft,fstft})';
% stft=griddata(tstft_rawM,fstft_rawM,ab,tstftM,fstftM);
% 
stft=stftRaw;%imgaussfilt(stftInterp,interpAmount_t/2);
% figure;subplot(2,1,1)
% imagesc(stftRaw)
% subplot(2,1,2)
% imagesc(stft)

%% Iterative Griffin and Lim algorithm


% istft=get_ispgm(lent,winInds,windowCenters,analysisWin,Fs,nIncs,S);
%
% figure;plot(real(istft)); hold on; %plot(imag(istft)); %plot(abs(istft));
% plot(SUT)
% legend('Real istft','real SUT')
winLenInterp=numel(fstft)*interpAmount_f;
winInterp=interp1(linspace(0,1,lent),win,linspace(0,1,winLenInterp));
winIndsInterp=(1:numel(fstft))-round(numel(fstft)/2);
overlapAmount=(sum(win))/winInc;%numel(winIndsInterp)/(windowCentersInterp(2)-windowCentersInterp(1)); % This is the "overlapamount" AFTER interpolation
analysisWin=winInterp/overlapAmount; % Analysis window for the inverse STFT

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


S0=sqrt(stft).*exp(1j*rand(size(stft))*2*pi);%.*(-1*(spgm<0));%.*exp(1j*rand(size(stft))*2*pi); % Seed spgm
S0Mag=abs(sqrt(stft)); % S0 has the correct amplitude, but not the correct phase
% % % % S0=sggm;
xt=get_ispgm(lent,winIndsInterp,windowCentersInterp,analysisWin,Fs,nIncsInterp,S0);
xt0=xt; % This is the first initial guess

maxIteration=200;
i=1;

% Convergence criterion
di=zeros(1,maxIteration);
diC=di;
diR=di;

h1=figure;
h1.Position=[55 117 990 861];%[-1528 112 821 865];
xlimsZoom=t(round(lent/2))+winLen*dt*[-2 2];

while i<maxIteration+1
    
    Si=get_spgm(nIncsInterp,winIndsInterp,windowCentersInterp,lent,winInterp,dt,winLenInterp,xt);         % Get STFT of present signal guess xt
    Sip1=S0Mag.*Si./abs(Si);%.*exp(1j*angle(Si));%.*Si./abs(Si);                           % Enforce magnitude along with calculated phase from Si
    Sip1(isnan(Sip1))=0;
    xt=get_ispgm(lent,winIndsInterp,windowCentersInterp,analysisWin,Fs,nIncsInterp,Sip1);
    % di(i)=sqrt(sum(sum(abs(abs(Si)-abs(Sip1)).^2))
    % di(i)=norm(abs(abs(Si)-abs(Sip1)).^2)/norm(abs(Si).^2)
%     di(i)=(norm(abs(Si)-abs(Sip1))).^2;%)/norm(abs(Si).^2)
    diC(i)=sqrt(    sum(sum( abs(  abs(S0) - abs(Si)  ).^2))  /    sum(sum(  abs(S0).^2 )) );%)/norm(abs(Si).^2)
%     diC2(i)=sqrt(    sum(sum( abs(  sqrt(abs(S0).^2) - sqrt(abs(Si).^2)  ).^2 ))  /    sum(sum(  abs(S0).^2 )) );%)/norm(abs(Si).^2)
%     diC(i)=sqrt(    norm( abs(  sqrt(abs(S0).^2) - sqrt(abs(Si).^2)  ).^2, 'fro')  /    norm(  abs(S0).^2 , 'fro') );%)/norm(abs(Si).^2)
    diR(i)=sqrt(sum(abs(SUT-xt).^2)/sum(abs(SUT).^2));%)/norm(abs(Si).^2)
    xc=max(abs(xcorr(SUT,xt,50)))/(sqrt(sum(abs(SUT).^2)*sum(abs(xt).^2)));
    diR(i)=(1-xc);%)/norm(abs(Si).^2)
    
    % plot(t,real(xt)); drawnow();
    i=i+1;
    
    if mod(i,5)==1
    updatePlot(h1,xlimsZoom,t,tstft,fstft,xt0,SUT,xt,stft,Si,diC,diR)
    end
end


    updatePlot(h1,xlimsZoom,t,tstft,fstft,xt0,SUT,xt,stft,Si,diC,diR)











function updatePlot(h1,xlimsZoom,t,tstft,fstft,xt0,SUT,xt,stft,Si,diC,diR)
figure(h1);
FS=16;
subplot(3,2,1)
imagesc(tstft,fstft,stft)
title('SUT STFT')
subplot(3,2,2)
imagesc(tstft,fstft,abs(Si).^2)
title('current Iteration STFT')
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





function istft=get_ispgm(lent,winInds,windowCenters,analysisWin,Fs,nIncs,S)
% Reconstruct temporal waveform from stft by the overlap and add technique

istft=zeros(1,lent); % Initialize variable

% figure;
for i=1:nIncs
    
    ift=nifft(S(:,i),Fs);
%     sutInds=mod(winInds+windowCenters(i),lent)+1;
    istft=istft+circshift(analysisWin,windowCenters(i)).*ift.';%,;
    
    
    
%     yyaxis left; hold off
% plot(abs(istft)); hold off
%     yyaxis right; 
%     plot(abs(ift)); hold on; 
%     plot(analysisWin/max(analysisWin)*max(abs(abs(nifft(S(:,i),Fs)))))
%     hold off
%     title(num2str(i))
%     drawnow()
end

% istft=circshift(istft,lent/2);
% istft=istft;
end






function spgm=get_spgm(nIncs,winInds,windowCenters,lent,win,dt,winLen,SUT)


spgm=zeros(lent,nIncs);

for i=1:nIncs
%     sutInds=mod(winInds+windowCenters(i),lent)+1;
    spgm(:,i)=nfft(SUT.*circshift(win,windowCenters(i)),dt);
end

end



















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
% imagesc(tstft,fstft,stft);
% ylim([-fmax fmax])
% xlabel('Time'); ylabel('Frequency');
% set(gca,'FontSize',FS)






