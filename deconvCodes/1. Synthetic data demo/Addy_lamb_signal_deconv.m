clear;clc

addpath('..\')
%% Input parameters

f = 300e3;                  % excitation frequency
ts = 0.1e-6;             	% time step (how often the signal is sampled)
fs = 1/ts;                	% sampling frequency
A = 10;                    	% excitation signal amplitude
N = 5;                     	% number of cycles for signal
len_measurment = 300e-6;   	% 300 microseconds - total time of signal recorded
n_events = 15;              % number of events        

%% Generating pulse

omega = 2*pi*f;
time = 0:ts:N/f;

V = A*(1-cos(omega*time/N)).*sin(omega*time).*(time<(N/f))/2;

%% plotting input signal

figure(1);clf
subplot(3,1,1)
plot(time*1e6,V);
axis tight;grid on;
xlim([0 len_measurment*1e6])
% xlabel('Time - [µs]');
ylabel('Amplitude - [V]')
title('Input signal')

%% time domain signal modelling using convolution

temp = (0:ts:len_measurment)';
time_events = zeros(length(temp),1);
occurances = randi(length(temp),n_events,1);
time_events(occurances) = 1;
time_events = time_events';

signal_conv = conv(V,time_events);
time_conv = (0:1:length(signal_conv)-1)'*ts;

%% plotting objective fn and convolved synthetic signal

subplot(3,1,2)
plot(temp*1e6,time_events);
axis tight;grid on;
% xlabel('Time - [µs]');
ylabel('Amplitude - [V]')
title('Synthetic Objective function')
ylim([0 1.2])

subplot(3,1,3)
plot(time_conv*1e6,signal_conv);
axis tight;grid on;
xlabel('Time - [µs]');
ylabel('Amplitude - [V]')
title('Convolved signal')

%% deconvolution using L1-MM


N    = length(V);     % Order of the filter = length of the exc signal
Hd = dfilt.dffir(V);  % to create the filter
[b,a] = tf(Hd);         % to get the transfer functions "a" and "b" from the filter

lam = 2;                                      % lam : regularization parameter
Nit = 50;                                     % Nit : number of iterations

tic
[x, cost] = deconvL1(signal_conv, lam, b, a, Nit);       % Run algorithm
runtime = toc;

conv_sig = conv(x,V,'same');
[~,delay] = max(envelope(V)); % time delay correction after convolution
conv_sig = vertcat(conv_sig(end-delay:end),conv_sig(1:end-delay-1));

%% plotting the results

figure(2);clf
subplot(2,1,1);hold on;
plot(temp*1e6,time_events,'b');
plot(time_conv*1e6,x,'r');
axis tight;grid on;box on
% xlabel('Time - [µs]');
ylabel('Amplitude - [V]')
title('Original objective function vs recovered')
ylim([-0.3 1.2])
legend('Original','Recovered')

subplot(2,1,2);hold on
plot(time_conv*1e6,signal_conv,'b');
plot(time_conv*1e6,conv_sig,'r');
axis tight;grid on;box on
xlabel('Time - [µs]');
ylabel('Amplitude - [V]')
title('Modeled signal vs convolved signal from recovered objective function')
legend('Original','Recovered')

%% plotting cost function

figure(3);clf
plot(1:Nit,cost)
xlabel('Iterations')
ylabel('Cost')
title('Cost vs Iterations')

%% 
rmpath('..\')