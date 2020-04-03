close all; clear;clc;

% Time domain construction of a multisine

% TASK 3.1. Time domain random phase multisine

N = 1000;   % Number of samples
K = 10;     % Number of excited frequencies

k = 1:N;
m = 1:K;

mn = m'.*k;

phi = 2*pi*rand([K,1]);     % phi is between 0 and 2pi

multisine = cos(2*pi*mn/N + phi);
multisine = sum(multisine,1);

% DFT of the multisine
multisine_fft = fft(multisine,[],2);
tresh = 10^(-10);
multisine_fft(abs(multisine_fft)<tresh) = 0;
multisine_fft = fftshift(multisine_fft);

fq = -N/2:N/2 -1;

figure; hold on;
plot(multisine.');
title('Multisine');

figure; hold on;
subplot(2,1,1);
stem(fq,abs(multisine_fft));
ylabel('|X|');
xlim([-50,50]);
subplot(2,1,2);
stem(fq,angle(multisine_fft));
ylabel('\angle X');
xlabel('Frequency (bin)');
xlim([-50,50]);

%% TASK 0.3.2. Frequency axis in Hz.

fs =100;

figure; hold on;
plot(k/fs,multisine.');
title('Multisine');
xlabel('time (s)')

figure; hold on;
subplot(2,1,1);
stem(fq*fs/N,abs(multisine_fft));
ylabel('|X|');
xlim([-10,10]);
subplot(2,1,2);
stem(fq*fs/N,angle(multisine_fft));
ylabel('\angle X');
xlabel('Frequency (Hz)');
xlim([-10,10]);

%% TASK 0.3.3. Excite specific frequency lines

omega = [4, 8, 12, 16, 20, 24];
omega = 2*pi*omega;
fs =200;
K = 6;
phi = 2*pi*rand([K,1]);     % phi is between 0 and 2pi
omegan = omega'.*k;

multisine = cos(omegan/fs + phi);
multisine = sum(multisine,1);

% DFT of the multisine
multisine_fft = fft(multisine,[],2);
tresh = 10^(-10);
multisine_fft(abs(multisine_fft)<tresh) = 0;
multisine_fft = fftshift(multisine_fft);

t = k/fs;

figure; hold on;
plot(t,multisine.');
title('Multisine');
xlabel('time (s)')

figure; hold on;
subplot(2,1,1);
stem(fq*fs/N,abs(multisine_fft));
ylabel('|X|');
xlim([-30,30]);
subplot(2,1,2);
stem(fq*fs/N,angle(multisine_fft));
ylabel('\angle X');
xlabel('Frequency (Hz)');
xlim([-30,30]);