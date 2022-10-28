clear; clc;
addpath('..\')

%% 10MHz incident signal

fs = 1000e6;
f = 10e6;
wavelength = fs/f;
sig = wavemaker(3.5, f, fs);
figure(1);clf; 
subplot(3,1,1);
plot(sig, 'LineWidth', 1); 
title(strcat('input signal, wavelength =',num2str(wavelength),' data pts'));
axis([0,3300,-2,2])

%% Sparse signal;
N = 3000;                    % N : length of signal
s = zeros(N,1);
k = [50,(50+wavelength*0.1), 500,(500+wavelength*0.6), 1200,(1200+wavelength*1.6), 2200,(2200+wavelength*4.6)];
s(k) = 1;
subplot(3,1,2);
plot(s, 'LineWidth', 1); 
title('distribution: objective function');
xlim([0 3300])

%% convoluted signal
y = conv(sig,s);
subplot(3,1,3);
plot(y, 'LineWidth', 1);hold on;
plot(abs(hilbert(y)), 'LineWidth', 1);
title('convoluted signal between input and districution');
xlim([0 3300])

%% deconvolution

N    = length(sig);     % Order of the filter = length of the exc signal
Hd = dfilt.dffir(sig);  % to create the filter
[b,a] = tf(Hd);         % to get the transfer functions "a" and "b" from the filter

lam = 1;                                      % lam : regularization parameter
Nit = 30;                                     % Nit : number of iterations

tic
[x, cost] = deconvL1(y, lam, b, a, Nit);       % Run algorithm
runtime = toc;

%% 

figure(2);clf
subplot(2,1,1);hold on
plot(s, 'LineWidth', 1); 
plot(x, 'LineWidth', 1)
title('Deconvolved objective function')
legend('Original','recovered')
box on

subplot(2,1,2);hold on
plot(y, 'LineWidth', 1)
plot(conv(sig,x), 'LineWidth', 1)
box on; axis tight
legend('Original signal','Deconvolved')
title('Signals')

%% 
rmpath('..\')

%% 

function x = wavemaker(nCycles, fc, fs)
% function to generate wave packet;
nSample = round(fs / fc * nCycles);
ts      = 1 / fs;
t_max   = ts * (nSample-1);
t       = 0: ts: t_max;
x       = sin( 2 * pi * fc .* t);
x = x.*hanning(nSample)';
end