addpath( '/Users/ben/Documents/MATLAB/library_repo' )


figDir='/Users/ben/Documents/MATLAB/timeFrequencyAnalysis/phaseRecoveryFinalGraphs/'

% fn='spgm_1  2  1  1  2_05h01m.fig';
fnt='temporalTrace_1  2  1  1  2_05h01m.fig';
fnZ='spgm_1  1  2  2_14h25m.fig';


%% Generate new data file
% % % % % open([figDir fnZ]);
%% oso
outputFigDir='/Users/ben/Desktop/PhD/manuscripts/Journals/TLS/phaseRecovery_TLS/figs/phaseRecovFigs/'

fn='spgm_1  1  1  1  1_14h33m.fig';zz=0; ylims=[-240 240];ylimPhase=[-1500 1850];
plotType='oso_'
open([figDir 'oso/' fn]);

% 
% % % % % zz
% fn='spgm_1  1  1  1_16h23m.fig';zz=1
% fn='spgm_1  1  1  1_17h16m.fig';zz=1; ylims=[-28.2 28.2];
% fn='spgm_1  1  1  1_18h36m.fig';zz=1; ylims=[-28.2 28.2]; ylimPhase=[-240 10];
% % fn='spgm_2  1  2  1_02h54m.fig';zz=1; ylims=[-28.2 28.2];
% % fn='spgm_1  1  2  2_02h29m.fig';zz=1; ylims=[-28.2 28.2];
% plotType='zz_'
% 
% h1=open([figDir 'zz/' fn]);





spgmOut=getimage(gcf);
h2=subplot(4,4,[3,4, 7,8])
delete(h2)
spgmIn=getimage(gcf);
h1=subplot(4,4,[1,2,5,6])
delete(h1)

[xdat,ydat]=getFigDat(gcf);

t=xdat{3};
yAbs=ydat{5};
yAngle=ydat{3};

% close all
%% process OSO data
% 
% load([figDir 'osoDat'])
% 
% spgmIn=osoInput;
% spgmOut=osoRecov; 
% 
% t=xdatOSO{3};
% yAbs=ydatOSO{5};
% yAngle=ydatOSO{3};
% ylims=[-226 226];


% figure;plot(x,yAbs)

%% process zigZagsignal
% 
% load([figDir 'zigZagDat'])

% spgmIn=zigzagInput;
% spgmOut=zigzagRecov; 
% 
% t=xdatZZ{3};
% yAbs=ydatZZ{5}; yAbs=yAbs/max(yAbs);
% yAngle=ydatZZ{3};
% ylims=[-28.2 28.2];
% 


% 
% spgmIn=spgmIn(:,1:2569);
% spgmOut=spgmOut(:,1:2569);
% t=t(1:10262);
% yAbs=yAbs(1:10262);
% yAngle=yAngle(1:10262);
% 

% Fix ymax regions

if zz
    iniInd=20; 
    finInd=0;
    t=t(iniInd:end-finInd);
    yAbs=yAbs(iniInd:end-finInd); yAngle=-yAngle(iniInd:end-finInd);
end
% get t-f vectors
dt=t(2)-t(1);
f=linspace(-1/(2*dt),1/(2*dt),numel(t))*1e-9;
df=f(2)-f(1);



y=yAbs.*exp(1j*yAngle);
yf=nfft(y);
yfAbs=abs(yf); yfAngle=unwrap(angle(yf));


dtSpec=(t(end)-t(1))/numel(spgmIn(1,:));
tSpec=(1:numel(spgmIn(1,:)))*dtSpec*1e9-dtSpec*1e9;
fSpec=linspace(-1/(2*dt),1/(2*dt),numel(spgmIn(:,1)))*1e-9;

% Fix vertical/frequency shift
if zz
    yfShiftGHz=2.611; % GHz
    yfShiftn=(-round(yfShiftGHz/df)-6);
    yfAbs=flip(circshift(yfAbs,-yfShiftn));
    yfAngle=flip(circshift(yfAngle,-yfShiftn)); yf=flip(circshift(yf,-yfShiftn));
    y=nifft(yf); yAbs=abs(y); yAngle=unwrap(angle(y));
    spgmIn=fliplr(circshift(spgmIn,yfShiftn));
    spgmOut=fliplr(circshift(spgmOut,yfShiftn));




zerPadLen=2^6;
shift=0;%-1.05e-9;
yForPadding=circshift(y,round(shift/dt));
ypad=[zeros(1,zerPadLen),yForPadding,zeros(1,zerPadLen)];
yfpad=nfft(ypad);dfPad=1/(numel(yfpad)*dt);
fpad=linspace(-1/(2*dt),1/(2*dt),numel(yfpad))*1e-9;

% Tried a peak search, didn't work very wekk
locs=1:numel(yfpad);
% [pks,locs,w,p]=findpeaks(abs(yfpad).^2,'MinPeakDistance',round(120e6/dfPad),...
%     'MinPeakHeight',0.2);
% nptPerLocs=3; locs=reshape((locs'+[-floor(nptPerLocs/2):floor(nptPerLocs/2)])',[numel(locs)*nptPerLocs 1]);
% locs=unique(locs);locs(locs<1)=[]; locs(locs>numel(yfpad))=[];

yAngleClean=unwrap(angle(yfpad(locs)));
yAngleClean=yAngleClean-(yAngleClean(round(numel(yAngleClean)/2)));

figure;plot(fpad,abs(yfpad)); hold on
plot(fpad(locs),abs(yfpad(locs)),'.');
yyaxis right ;
% plot(unwrap(angle(yfpad)))
% hold on
plot(fpad(locs),yAngleClean)
ylim([-100 200])
xlim(ylims)

end


mtIn=sum(spgmIn);mtIn=mtIn/max(mtIn);
mtOut=sum(spgmOut);mtOut=mtOut/max(mtOut);
mfIn=sum(spgmIn,2);mfIn=mfIn/max(mfIn);
mfOut=sum(spgmOut,2);mfOut=mfOut/max(mfOut);



figure;
imagesc(tSpec,fSpec,spgmIn/(max(max(spgmIn))))
colormap('jet'); ylim(ylims)
title('deconvolved spgm'); xlabel('Time (ns)'); ylabel('Frequency (GHz)');
save_figsvgeps([outputFigDir plotType 'spgmIn'])

figure
imagesc(tSpec,fSpec,spgmOut/(max(max(spgmOut))))
colormap('jet');ylim(ylims)
title('reconstructed spgm'); xlabel('Time (ns)'); ylabel('Frequency (GHz)');
save_figsvgeps([outputFigDir plotType 'spgmOut'])

figure
% subplot(2,2,3)
plot(t*1e9,yAbs.^2/max(yAbs.^2));
hold on; 
% plot(tSpec,mtIn)
plot(tSpec,mtOut)
xlim([tSpec(1) tSpec(end)])
legend('Abs','mtIn','mtOut')
ylabel('Intensity (n.u.)');
yyaxis right;
% plot(t*1e9,yAngle);
if zz
    phase=yAngle;
plot(t*1e9,phase-max(phase))
else
    plot(t*1e9,yAngle-yAngle(round(numel(yAngle)/2))); ylim(ylimPhase)
end
ylim(ylimPhase)
ylabel('Phase (rad)');
title('reconstructed temporal trace'); xlabel('Time (ns)'); 
save_figsvgeps([outputFigDir plotType 'time'])





figure
% subplot(2,2,4)
plot(f,yfAbs.^2/max(yfAbs.^2)); 
hold on; 
% plot(fSpec,mfIn)
plot(fSpec,mfOut)
legend('Abs','mfIn','mfOut')
% 
% hold on;
% plot(fSpec,sum(spgmOut,2)/max(sum(spgmOut,2)))

yyaxis right
if zz
    phase=unwrap(angle(yfpad(locs)));
plot(fpad(locs),phase-phase(round(numel(phase)/2))-60);%max(phase))
else
    plot(f,yfAngle-yfAngle(round(numel(yfAngle)/2))); ylim(ylimPhase)
end
xlim(ylims);ylim(ylimPhase)
ylabel('Phase (rad)');
title('reconstructed spectral trace'); xlabel('Frequency (GHz)'); 
save_figsvgeps([outputFigDir plotType 'freq'])
