% https://eeweb.engineering.nyu.edu/iselesni/lecture_notes/sparse_deconv/index.html

clear;clc
addpath('..\')

load exc_375.mat;    % loading exitation signal (input signal or the point spread function (PSF))
load s0a0bscan.mat;  % loading the recorded ultrasoinc signals
load deconv_bscan_ivan;     % load already calculate reflectvity for original data

% load deconv_bscan_10SNR_ivan.mat  % load already calculate reflectvity for noisy data
% Bscan = awgn(Bscan,10);   % this is how i added noise to my data

time = exc(:,1);    % time vector for x-axis
exc(:,1)=[];        % removing time vector
exc(268:end) = [];  % cropping the zero padding of the input signal

%% Create data

fs = 1/(time(2)*1e-6);  % Sampling Frequency from the time vector

N    = length(exc);     % Order of the filter = length of the exc signal
Hd = dfilt.dffir(exc);  % to create the filter
[b,a] = tf(Hd);         % to get the transfer functions "a" and "b" from the filter
h = filter(b, a, [1 zeros(1,N)]);       % Calculate impulse response. Just to check if you get the same PSF

figure(1);clf
plot(time(1:N+1),h)
box on;axis tight;grid on
title('Excitation function')
xlabel('Time - [µs]')
ylabel('Amplitude')

%% Plot the B-Scan signal

figure(2);clf
pcolor(time,(1:200),abs(Bscan'));
shading interp
colormap hot
colorbar
xlabel('Time - [µs]')
ylabel('Distance - [mm]')
title('Bscan')

%%%%----------------------------------
% Uncomment from next line to next dash line at line 94 
% for demo of one Lamb wave signal and comment rest of the code till end
% %% 
% 
% lam = 2;                                      % lam : regularization parameter
% Nit = 20;                                     % Nit : number of iterations
% 
% i = 100;                                      % which signal to load from Bscan
% y = Bscan(:,i);
% 
% if exist('x','var')
%     clear x
% end
% 
% [x, cost] = deconvL1(y, lam, b, a, Nit);       % Run algorithm
% 
% rmpath('..\')
% %% 
% 
% figure(3);clf
% subplot(2,1,1); hold on;
% plot(time,y);
% plot(time,conv(x,h,'same'));
% box on;axis tight;grid on
% title('A-Scan')
% xlabel('Time - [µs]')
% ylabel('Amplitude')
% legend('Output signal y(n)','Deconvolved signal')
% ylim([min(min(Bscan)) max(max(Bscan))])
% 
% temp = x;
% subplot(2,1,2)
% plot(time,temp);
% box on;axis tight;grid on
% title('Reflectivity sequence using l1- MM norm')
% xlabel('Time - [µs]')
% ylabel('Amplitude')
% ylim([min(min(x)) max(max(x))])
% 
% sgtitle(['Signal at D = ' num2str(i) ' mm'])
% 
% %% 
% 
% figure(5)
% clf
% plot(cost,'linewidth',2)
% title('Cost function history')
% xlabel('Iterations')
% box on;grid on
% xlim([0 Nit])
%%%%----------------------------------

%% Sparse deconvolution for all the B-scans

lam = 2;                                      % lam : regularization parameter
Nit = 20;                                     % Nit : number of iterations

skip = exist('x','var');

if ~skip
    f = waitbar(0,'calculating');
    x = zeros(size(Bscan));
    cost = zeros(size(Bscan,2),Nit);
    for i = 1:size(Bscan,2)
        
        waitbar(i/size(Bscan,2),f,[num2str(i) ' out of ' num2str(size(Bscan,2))]);
        y = Bscan(:,i);
        
        [x(:,i), cost(:,i)] = deconvL1(y, lam, b, a, Nit);       % Run algorithm
        
    end
    close(f)
end

%% 

i = 1;
figure(3);clf
subplot(2,1,1); hold on;
gca1 = plot(time,Bscan(:,i));
gca2 = plot(time,conv(x(:,i),h,'same'));
box on;axis tight;grid on
title('A-Scan')
xlabel('Time - [µs]')
ylabel('Amplitude')
legend('Output signal y(n)','Deconvolved signal')
ylim([min(min(Bscan)) max(max(Bscan))])

temp = x(:,i);
subplot(2,1,2)
gca3 = plot(time,temp);
box on;axis tight;grid on
title('Reflectivity sequence using l1- MM norm')
xlabel('Time - [µs]')
ylabel('Amplitude')
ylim([min(min(x)) max(max(x))])

sgtitle(['Signal at D = ' num2str(i) ' mm'])

for i = 1:size(Bscan,2)
    
    set(gca1,'YData',Bscan(:,i))
    set(gca2,'YData',conv(x(:,i),h,'same'))
    set(gca3,'YData',x(:,i))
    sgtitle(['Signal at D = ' num2str(i) ' mm'])
    drawnow
    
end


%% 

figure(4);clf
pcolor(time,(1:200),abs(x'));
shading interp
colormap hot
colorbar
xlabel('Time - [µs]')
ylabel('Distance - [mm]')
title('Deconvolved Bscan')

rmpath('..\')