% https://dsp.stackexchange.com/questions/71822/deconvolution-of-a-1d-time-domain-wave-signal-convolved-with-series-of-rect-sign
% https://github.com/RoyiAvital/StackExchangeCodes/tree/master/SignalProcessing/Q71822

clear;clc
addpath('..\')
load exc_375.mat;    % loading exitation signal (input signal or the point spread function (PSF))
load s0a0bscan.mat;  % loading the recorded ultrasoinc signals
load deconvbscan;

time = exc(:,1);    % time vector for x-axis
exc(:,1)=[];        % removing time vector
vH = exc;           % just storing in a differnent variable for matching Royi's code
vH(268:end,:)=[];   % cropping the zero padding of the input signal

%% Simulation Parameters

convShape       = 1;
paramLambda     = 1;
numIterations   = 1000;

%% This and the next section is for deconvolving one signal and seeing the results

% i = 105; % index of the signal to be loaded 
% vY = Bscan(:,i);    % loading the signal
% numElementsX = size(vY, 1) - size(vH, 1) + 1;   % length of signal after deconvolution
% mH = CreateConvMtx1D(vH, numElementsX, convShape);  % sparse matrix for input signal
% mD = CreateConvMtx1D([1, 0], numElementsX, convShape); % sparse matrix - but I dont know what is this for
% vXHat   = SolveLsTvAdmm([], full(mH), vY, full(mD), paramLambda, numIterations); %L1 ADMM solver
% 
% %% 
% 
% temp = vXHat;
% % temp(abs(temp)>2)=0;        % thresholding to remove big spikes
% 
% figure(1);clf
% subplot(2,1,1);hold on
% plot(time,vY)
% plot(time,conv(temp,vH))
% box on;grid on
% xlim([0 time(end)])
% xlabel('Time - [µs]')
% ylabel('Amplitude')
% legend('Oridignal signal','Convolved from reflectivity sequence')
% 
% 
% subplot(2,1,2)
% plot(time(1:length(temp)),temp)
% box on;grid on
% xlim([0 time(end)])
% title('Reflectivity sequence')
% 
% %% Plot the B-Scan signal
% 
% figure(2);clf
% pcolor(time,(1:200),abs(Bscan'));
% shading interp
% colormap hot
% colorbar
% xlabel('Time - [µs]')
% ylabel('Distance - [mm]')
% title('Bscan')

%% comment previous two sections and Uncomment this section for deconvolution of all A-scans

skip = exist('vXHat','var');
if ~skip
    f = waitbar(0,'calculating');

    for i = 1:size(Bscan,2)
        vY = Bscan(:,i);
        waitbar(i/size(Bscan,2),f,[num2str(i) ' out of ' num2str(size(Bscan,2))]);
        numElementsX = size(vY, 1) - size(vH, 1) + 1;
        mH = CreateConvMtx1D(vH, numElementsX, convShape);
        mD = CreateConvMtx1D([1, -1], numElementsX, convShape);
        vXHat(:,i)   = SolveLsTvAdmm([], full(mH), vY, full(mD), paramLambda, numIterations);
    end

    close(f)
end

%% ploting

temp = vXHat;
% temp(abs(temp)>0.8)=0;        % thresholding to remove big spikes

figure(100);clf
pcolor(time(1:length(temp)),1:200,abs(temp'))
xlabel('Time - [µs]')
ylabel('Distance - [mm]')
shading interp; axis tight
colormap(flipud(gray))
colorbar
caxis([0.001 0.009])
title('Reflectivity sequence B-Scan')


rmpath('..\')