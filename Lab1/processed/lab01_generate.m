clear
close all

Fs = 8e3;
FsdivN = 1;
N = Fs/FsdivN;
des_rms = 1e-3; %1mV

K = 500;
t = (0:N-1)/Fs;
f = (0:N-1)*Fs/N;

A = ones(1,K);
k = 0:K-1;
phi = k.*(k+1)*pi/K; %Schroeder
phi = zeros(1,K); %zero;
phi = rand(1,K)*2*pi; %random
X = zeros(1,N);
X(2:K+1) = A.*exp(1j*phi);
x = N*real(ifft(X));
x = x*des_rms/rms(x);
X = fft(x);
CF = max(abs(x))/rms(x);

figure
subplot(2,1,1)
plot(t,x)
xlabel('Time (s)')
ylabel('Voltage (V)')
title(['x(t), RMS = ' num2str(rms(x)*1000) 'mV, CF = ' num2str(CF)])
subplot(2,1,2)
plot(f/1000,abs(X),'.--')
xlabel('Frequency (kHz)')
title('|DFT(x)|')


save('Group1_Input1.mat','u_sch_500')

