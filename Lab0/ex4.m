close all;clear;clc;

% 0.4 Frequency domain construction of a multisine

N = 1000;
m = 30;
phi = 2*pi*rand([m,1]);

X = zeros([N,1]);

X(1:m) = exp(1i*phi);

x = ifftshift(N*real(ifft(X)));

figure; hold on;
plot(x);


figure; hold on;
subplot(2,1,1);
stem(abs(X));
ylabel('|X|');
xlim([-50,50]);
subplot(2,1,2);
stem(angle(X));
ylabel('\angle X');
xlabel('Frequency (bin)');
xlim([-50,50]);


% TASK 0.4.3. Specified excited frequency band and frequency resolution

phi = 2*pi*rand([m,1]);
r = 3;

X = zeros([N,1]);
X(5*3:15*3-1) = exp(1i*phi);

figure; hold on;
subplot(2,1,1);
stem(abs(X));
ylabel('|X|');
xlim([-50,50]);
subplot(2,1,2);
stem(angle(X));
ylabel('\angle X');
xlabel('Frequency (bin)');
xlim([-50,50]);