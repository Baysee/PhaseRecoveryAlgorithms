% % % Vanilla Frog Algorithm


% Spectrogram function
addpath( '/Users/ben/Documents/MATLAB/library' )
%% Time frequency vectors definition

lent=2^19;                      % Signal length
tWind=50e-9;                   % Time window span

% t=linspace(0,tWind,lent);
t=linspace(-tWind/2,tWind/2,lent);
dt=t(2)-t(1);Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
fGHz=f*10^-9;tps=t*10^12;%GHz
fG=1e-9; tn=1e9;
scale=1;

spectrogramTypes={'T-TAI','Time-lens'};
processUsingType=2; % Make =1 for T-TAI, =2 for Time lens
spectrgramNow=spectrogramTypes{processUsingType};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% temporal phase signal definiton %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These settings are used for both the T-TAI and the time-lens spectrogram
trTry=0.95e-9;%9.0761e-9;        % Time-lens apperture
eta=20;                % number of analysis points
fmax=eta/trTry;         % full bandwidth that can be processed with the spectrogram

ntr=round(trTry/dt);    % tr is adjusted to be an integer amount of points
tr=ntr*dt;
nSamps=ceil(lent/ntr);
tSingleSamp=dt*(1:ntr);
lensInds=0:ntr:lent;


switch spectrgramNow
    
    case 'Time-lens'
        %% Time Lens Approach
        
        C_tl=2*pi*eta/(tr^2);
        singleSamp=C_tl/2*(tSingleSamp-tSingleSamp(round(ntr/2))).^2;
        sampSigRaw=repmat(singleSamp,1,nSamps);
        sampSigRawLent=sampSigRaw(1:lent);
        sampSig=real(filtSG_tf(sampSigRawLent,t,f,round(800e9/df),5,0));
        sampSig=sampSig(1:lent);
        sampFreq=1/tr
        tr*C_tl/(2*pi)
        
        %%%%%%%%%%%%%%%%
        %% Dispersion %%
        %%%%%%%%%%%%%%%%
        phi=1/C_tl;
        
    case 'T-TAI'
        
        %% Talbot TAI - Uncomment below for TAI and comment "Time Lens"
        
        m=eta; % number of analysis points is equivalent to "amplification factor"
        ts=tr/m;
        AWG_nuq=1/ts;
        p=1;
        s=generateSparameter(p,m);
        GV=wrapTo2Pi(s/m*pi*((0:m-1)).^2);
        nSampsPer_ts=ceil(ntr/m);
        GVtry=repelem(GV,1,nSampsPer_ts);
        t_samples=(1:numel(GVtry))/numel(GVtry)*tr;
        singleSamp=interp1(t_samples,GVtry,tSingleSamp);
        allSamps=repmat(singleSamp,1,nSamps);
        allSamps=allSamps(1:lent);
        sampSig=real(filtSG_tf(allSamps,t,f,round(200e9/df),10,1));
        
        %%%%%%%%%%%%%%%%
        %% Dispersion %%
        %%%%%%%%%%%%%%%%
        
        phi=p*m*(tr/m)^2/(2*pi);% this is equal to p*m*ts^2/(2*pi)
        
end


% Leave below uncommented for all cases.

phi2perKm=   2.1823e-23;
totalDispersion_equivalent_SMF_km=phi/phi2perKm;
phaseGVD=phi/2*(2*pi*f).^2;



%%%%%%%%%%%%%%%%%%
%% Generate SUT %%
%%%%%%%%%%%%%%%%%%
SUT=ones(1,lent);



% Linearly chirped signal
fini=-fmax/2;0; ffin=fmax/2;
c=(ffin-fini)/tWind;
% SUT=sin(2*pi*(c/2*t.^2+fini*t));

sutWind=tWind/16; sutCenter=(t(end)+t(1))/2;
SUT=ones(1,lent).*exp(1j*2*pi*(c/2*t.^2+fini*t));
SUT=superGauss(0,sutWind,10,t,sutCenter);

% Help setup frequency range by carving slightly bigger than sut, then
% multiply by gaussain frequency chirp
% superGauss(0,sutWind*1.5,30,t,t(end)/2).*

% instantaneous frequency parameters
fEdge=fmax/2/1.5
fMid=-fmax/2/1.6
centerMid=sutCenter+tWind/50;
instF_sigma=sutWind/1.8;
SUTinstFIni=superGauss(0,instF_sigma,1,t,centerMid);
riseInd=find(SUT>0.5*max(SUT),1); fallInd=find(SUT>0.5*max(SUT),1,'last'); 
% rangeY=max(SUTinstFIni)-min([SUTinstFIni(riseInd), SUTinstFIni(fallInd)]);
% SUTinstF=SUTinstFIni*(abs(fEdge-fMid))/rangeY+min(fEdge,fMid);
extraBitFactor=20;
SUTregion=riseInd-round(sutWind/dt/extraBitFactor)+1:-1+fallInd+round(sutWind/dt/extraBitFactor);
extraBitFactorPlus=4;US=20;
SUTregP_US=riseInd-round(sutWind/dt/extraBitFactorPlus)+1:US:-1+fallInd+round(sutWind/dt/extraBitFactorPlus);

SUTinstFIni([1:SUTregion(1)-1,SUTregion(end)+1:end])=0;
SUTinstFIni(SUTregion)=SUTinstFIni(SUTregion)-min(SUTinstFIni(SUTregion)); SUTinstFIni=SUTinstFIni/max(SUTinstFIni);
SUTinstF=SUTinstFIni*(abs(fEdge-fMid))+min(fEdge,fMid);
SUTphase=cumtrapz(t,SUTinstF);

SUT=SUT.*exp(1j*2*pi*SUTphase);



figure;
subplot(3,1,1)
plot(t,abs(SUT).^2); hold on ; plot(t,SUTinstFIni); 
yyaxis right
plot(t,unwrap(angle(SUT)))
xlabel('Time')
legend('sut power','sutinstFini','anle SUT')
subplot(3,1,2)
plot(SUTinstF); 
yyaxis right; plot(SUTphase)
legend('SUTinstF','SUTphase')
xlabel('Time')

subplot(3,1,3)
plot(f,abs(nfft(SUT)))
xlabel('Freq')
xlim([-fmax/2 fmax/2])


%%%%%%%%%%%%%%%%
%% Processing %%
%%%%%%%%%%%%%%%%


sampSUT=exp(1j*sampSig).*SUT   ;
sampSUT_f=nfft(sampSUT,dt,scale);

sampSUTdisp_f=(sampSUT_f).*exp(1j*phaseGVD);
sampSUTdisp=nifft(sampSUTdisp_f,Fs,scale);


%% Filter Spectrogram Signal
sampSUTdisp=(filtSG_tf(sampSUTdisp,t,f,round((10/tr)/df),5,0));

%%% Reshape output signal to get 2D representation
spgm=reshape(circshift(sampSUTdisp(1:lensInds(end)),0),ntr,(numel(lensInds)-1));
f_spec=(tSingleSamp-max(tSingleSamp)/2)/phi/(2*pi);
t_spec=linspace(-tWind/2,tWind/2,numel(lensInds)-1);



%% Plotting section

ylims=([-0.1 1.1]);
ylimsNegPos=([-1.1 1.1]);
numtrWindLims=20;
tlims=([-20 20]);%'auto';%(t(end)/2+[-numtrWindLims*tr numtrWindLims*tr])*1e9;
% tlims=[t(1) t(end)];

tps=1e12; tn=1e9; tus=1e6;
lbps='Time (ps)'; lbns='Time (ns)'; lbus='Time (us)';
fG=1e-9; lbG='Frequency (GHz)';
figure('Renderer', 'painters', 'Position', [50 50 1200 700])




subplot(4,4,1:4)
plot(t*tn,abs(SUT).^2/(max(abs(SUT).^2)))
hold on
plot(t*tn,abs(sampSUTdisp).^2/(max(abs(sampSUTdisp).^2)))
ylabel('Intensity')
yyaxis right
plot(t*tn,sampSig)
%  xlim([t(1) t(end)]*tns)
xlim(tlims);%ylim(ylims);
xlabel(lbns)
legend('sampled SUT','Dispersed wvf','Time Lens')
ylabel('Phase (rad)')



subplot(4,4,5:8)
plot(t*tn,abs(SUT));%/(max(abs(SUT))))
hold on
plot(t*tn,angle(SUT));%/(max(angle(SUT))))

plot(t*tn,abs(sampSUTdisp).^2/(max(abs(sampSUTdisp).^2)))
yyaxis right
plot(t*tn,sampSig)

ylabel('Intensity')
% yyaxis right
% plot(t*tns,sampSig)
xlim(tlims);%ylim(ylimsNegPos);
xlabel(lbns)
legend('abs','angle','spectrogram')
ylabel('Phase (rad)')

subplot(4,4,9:12)
% plot(t*tns,real(SUT)/(max(real(SUT))))
% hold on
% plot(t*tns,imag(SUT)/(max(imag(SUT))))

plot(t*tn,abs(sampSUTdisp).^2/(max(abs(sampSUTdisp).^2)))
ylabel('Intensity')
% yyaxis right
% plot(t*tns,sampSig)
xlim(tlims);%ylim(ylimsNegPos);
xlabel(lbns)
legend('real','imag','spectrogram')
ylabel('Phase (rad)')

subplot(4,4,13:16)

imagesc(t_spec*tn,f_spec*fG,(abs(spgm).^2));
hcb=colorbar();
xlim(tlims)

ylabel(hcb,'Power (dB)')
xlabel('Time (ns)')







%% Detailed color plots


% samp SUT instantaneous frequency

% cmap=[(logspace(0,1,nCmapElem))'/10,zeros(nCmapElem,1)+0.2,(logspace(1,0,nCmapElem))'/10];

sutrp=abs(SUT(SUTregP_US)).^2;
sutIF=gradient(unwrap(angle(SUT)),dt)/(2*pi);
sutIF=sutIF(SUTregP_US);
sutModIF=gradient(unwrap(angle(sampSUT)),dt)/(2*pi);
sutModIF=sutModIF(SUTregP_US);
zerosSRP=zeros(size(sutIF))+min(sutrp);
tIF=t((SUTregP_US));


% Output signal instantaneous frequency
% pksInds=abs(sampSUTdisp).^2>max(abs(sampSUTdisp).^2)*0.5;
[pks,locs]=findpeaks(abs(sampSUTdisp),'MinPeakDistance',ntr*0.8,'MinPeakHeight',0.2*max(abs(sampSUTdisp)));
ifPeakRange=round(ntr/eta/4);
locs2d=locs'+[-round(ifPeakRange/2):round(ifPeakRange/2)];
peakInds=reshape(locs2d',[1,numel(locs2d)]);
 
outputInstF=1/(2*pi)*gradient(unwrap(angle(sampSUTdisp)),dt);
outputPeakInstF=outputInstF(peakInds);
peakIndsinterp=[SUTregP_US(1),peakInds,SUTregP_US(end)];
outputPeakInstFinterp=[outputPeakInstF(1),outputPeakInstF,outputPeakInstF(end)];
outputInstF_interp=interp1(peakIndsinterp,outputPeakInstFinterp,SUTregP_US);

outputIntensityRP=abs(sampSUTdisp(SUTregP_US)).^2;
outputIF=outputInstF_interp;

% plot all

nCmapElem=200;
% cmap=[linspace(0.2,0.85,nCmapElem)',zeros(nCmapElem,1)+0.1,linspace(0.85,0.2,nCmapElem)'];
% cmap=[linspace(0.1,0.85,nCmapElem)',zeros(nCmapElem,1)+0.05,linspace(0.9,0.1,nCmapElem)'];
% cmap=jet(1000);
 cmap=[linspace(0.02,0.98,nCmapElem)',zeros(nCmapElem,1)+0.05,linspace(0.98,0.02,nCmapElem)'];
cmp=cmap;
% cmp = getPyPlot_cMap('CMRmap',1400);%,keepAlpha,pyCmd);
% cmp=cmp*0.9
figure;

xti='';
yti='';
zti='';

% SUT MOD 

figure
% s=surf([tIF;tIF],[zerosSRP;sutrp],[zerosSRP;zerosSRP],[sutIF;sutIF]);
s=surf([tIF;tIF],[zerosSRP;zerosSRP],[zerosSRP;sutrp],[sutModIF;sutModIF]);
s.EdgeColor='none'
colormap(cmp)
% az=-39.9; el=21.2;
az=0; el=0;
view(az,el)
caxis2=caxis*0.9;
caxis(caxis2)
xticks(xti); yticks(yti); zticks(zti);

% SUT

figure

s=surf([tIF;tIF],[zerosSRP;zerosSRP],[zerosSRP;sutrp],[sutIF;sutIF]);
s.EdgeColor='none'

colormap(cmp)
caxis(caxis2);
view(az,el)
xticks(xti); yticks(yti); zticks(zti);



% output 

figure;
s=surf([tIF;tIF],[zerosSRP;zerosSRP],[zerosSRP;outputIntensityRP],[outputIF;outputIF]);
s.EdgeColor='none'
colormap(cmp)
 view(az,el)
caxis(caxis2);
ylim([-1 3])
xticks(xti); yticks(yti); zticks(zti);


% Lens
figure
verticalPart=[(sampSig(SUTregP_US)-(max(sampSig(SUTregP_US))-min(sampSig(SUTregP_US)))*0.1);1.3*sampSig(SUTregP_US)];
s=surf([tIF;tIF],...
    [zerosSRP;zerosSRP],...
    verticalPart,...
    [gradient(sampSig(SUTregP_US));gradient(sampSig(SUTregP_US))]);
colormap(cmp)
s.EdgeColor='none'
caxisLens=caxis;
caxis(caxisLens*1.2)
 view(az,el)

zlim([min(min(verticalPart)) max(max(verticalPart))])
xticks(xti); yticks(yti); zticks(zti);



% Spectrogram plot


% imagesc(t_spec*tns,f_spec*fG,(abs(spgm).^2));

figure;
[t2d,f2d]=meshgrid(t_spec*tn,f_spec*fG);
zVals=zeros(size(t2d));
figure
s=surf([t2d],...
    [f2d],...
    zVals,...
    (abs(spgm).^2));
% zlim([-0.1 0.1])
xlim(1e9*[tIF(1) tIF(end)]+tr*1e9)
ylim([f_spec(1)*fG f_spec(end)*fG])
s.EdgeColor='none'
view(0,90)
xticks(xti); yticks(yti); zticks(zti);



figure;
imagesc(t_spec*tn,f_spec*fG, (abs(spgm).^2))
xlim(1e9*[tIF(1) tIF(end)]+tr*1e9)
ylim([f_spec(1)*fG f_spec(end)*fG])
