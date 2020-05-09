clear; close all;clc;
%%
% 
%  Measuring transfer functions
%  
% 

% Generation of multisine signals.

t = 1:4096;     % Time
K_max = 100;
K = 1:K_max;      % Number of excited spectral lines

% Constant phase multisine
phi = 0;
multisine = sum(sin(2*pi*t'*K/length(t) + phi),2);

% Random phase mutlisine
phi = 2*pi*rand([1,K_max]);
multisine_rand = sum(sin(2*pi*t'*K/length(t) + phi),2);

% Schroeder phase
phi = K.*(K+1)*pi/length(K);
multisine_schro = sum(sin(2*pi*t'*K/length(t) + phi),2);

% Crest factor

% Crest factor - constant phase multisine
crestFact = max(multisine)/rms(multisine); % rms requires the signal
%Toolbox
% crestFact = max(multisine)*sqrt(length(multisine))/norm(multisine);

% Crest factor - random phase mutlisine
crestFact_rand = max(multisine_rand)/rms(multisine_rand);
% crestFact_rand = max(multisine_rand)*sqrt(length(multisine_rand))/norm(multisine_rand);

% Crest factor - Schroeder phase mutlisine
crestFact_schro = max(multisine_schro)/rms(multisine_schro);
% crestFact_schro = max(multisine_schro)*sqrt(length(multisine_schro))/norm(multisine_schro);

disp('Crest factors: ');
disp(join(['Constant phase multisine:  ',num2str(crestFact)]));
disp(join(['Random phase multisine:    ',num2str(crestFact_rand)]));
disp(join(['Schroeder phase multisine: ',num2str(crestFact_schro)]));

% FFTs

% FFT - constant phase multisine
fft_const = fftshift(fft(multisine));

% FFT - random phase multisine
fft_rand = fftshift(fft(multisine_rand));

% FFT - Schroeder phase multisine
fft_schro = fftshift(fft(multisine_schro));

% % Plots
% fq = -length(t)/2:length(t)/2 -1;
% 
% figure('Position',[500 500 800 420]);hold on;
% subplot(1,2,1)
% plot(t',multisine',t',multisine_rand,t',multisine_schro);
% legend('Constant phase spectrum','random phase spectrum','Schroeder phase spectrum');
% 
% subplot(1,2,2)
% plot(fq,abs(fft_const),fq,abs(fft_rand),fq,abs(fft_schro));
% xlim([-300 300]);
% legend('Constant phase spectrum','random phase spectrum','Schroeder phase spectrum');
% 
% figure;hold on;
% subplot(3,2,1);
% plot(t',multisine');
% xlabel('Time samples');
% 
% subplot(3,2,3);
% plot(t',multisine_rand');
% xlabel('Time samples');
% 
% subplot(3,2,5);
% plot(t',multisine_schro');
% xlabel('Time samples');
% 
% subplot(3,2,2);
% plot(fq,abs(fft_const));
% xlabel('Frequency');
% 
% subplot(3,2,4);
% plot(fq,abs(fft_rand));
% xlabel('Frequency');
% 
% subplot(3,2,6);
% plot(fq,abs(fft_schro));
% xlabel('Frequency');

%%
%
% DUT measurements
%
%

% DUT parameters

bandwidth = 500;
fs = 8e3;          % Sampling frequency

% Generation of the multisine 

%%
% $\omega_{0} = \frac{2\pi f_{s}}{N}$

N = fs;
t = 1:N;
K_max = 500;
K = 1:K_max;      % Number of excited spectral lines

% Schroeder phase

% Frequency domain construction
phi = K.*(K+1)*pi/length(K);
X = zeros(N,1);
X(2:K_max+1) = exp(1j*phi);
x = N*real(ifft(X));
% u_sch_500 = x*0.1/rms(x);
u_sch_500 = x*sqrt(length(x))/norm(x);

%FFT
fft_u_sch_500 = fftshift(fft(u_sch_500));

% Plots
% figure('Position',[500 500 800 420]);hold on;
% subplot(1,2,1)
% plot(t',u_sch_500);
% legend('Schroeder phase spectrum');
% 
fq = -length(t)/2:length(t)/2 -1;
% subplot(1,2,2)
% plot(fq,abs(fft_u_sch_500));
% %plot(fq,abs(X));
% xlim([-700 700]);
% legend('Schroeder phase spectrum');
% 
% % Save the input signal
% save('Group9_Input1.mat','u_sch_500');

%%
% Excitation signals
%

% Periodic

% Constant phase multisine
phi = 2*pi*0;
multisine_const = sum(sin(2*pi*t'*K/length(t) + phi),2);
multisine_const = 0.1*multisine_const/rms(multisine_const);
% multisine_const = 0.1*multisine_const*sqrt(length(multisine_const))/norm(multisine_const);

% Random phase mutlisine
phi = 2*pi*rand([1,K_max]);
multisine_rand = sum(sin(2*pi*t'*K/length(t) + phi),2);
multisine_rand = 0.1*multisine_const/rms(multisine_rand);
% multisine_rand = 0.1*multisine_const*sqrt(length(multisine_rand))/norm(multisine_rand);

% Periodic noise signal
d = 1;
period = N/d;
noise_per = randn([period,1]);
noise_per = repmat(noise_per,d);
noise_per = 0.1*noise_per/rms(noise_per);
% noise_per = 0.1*noise_per*sqrt(length(noise_per))/norm(noise_per);

% Aperiodic

noise_aper = randn([N,1]);
noise_aper = 0.1*noise_aper/rms(noise_aper);
% noise_aper = 0.1*noise_aper*sqrt(length(noise_aper))/norm(noise_aper);

% Multiplyin the noise by the Hanning window
[Hann] = hanning(N,'periodic');
noise_aper_h = noise_aper.*Hann;

figure('Position',[500 500 800 420]);hold on;
subplot(1,2,1);
plot(noise_aper);
subplot(1,2,2);
plot(t,noise_aper_h,t,Hann);


%%
addpath('./processed/');
load('Group09_Output1.Mat')
data = load('Group09_Output1.Mat');

figure;
subplot(10,2,1)
plot(Su);
title('Su')
subplot(10,2,2)
plot(Sy);
title('Sy')
for i = 1:9
    subplot(10,2,i*2+1)
    u = sprintf('Su%d',i);
    plot(data.(u));
    title(join(['Su',num2str(i)]));
    subplot(10,2,i*2+2)
    y = sprintf('Sy%d',i);
    plot(data.(y));
    title(join(['Sy',num2str(i)]));
end

figure;
subplot(10,2,1)
plot(fq,abs(fftshift(fft(Su))));
xlim([-510 510]);
title('Su')
subplot(10,2,2)
plot(fq,abs(fftshift(fft(Sy))));
xlim([-510 510]);
title('Sy')
for i = 1:9
    subplot(10,2,i*2+1)
    u = sprintf('Su%d',i);
    plot(fq,abs(fftshift(fft(data.(u)))));
    xlim([-510 510]);
    title(join(['Su',num2str(i)]));
    subplot(10,2,i*2+2)
    y = sprintf('Sy%d',i);
    plot(fq,abs(fftshift(fft(data.(y)))));
    xlim([-510 510]);
    title(join(['Sy',num2str(i)]));
end


frf9 = FRF(Su9,Sy9,5);
fq_small = -500:500-1;

figure('Name','output1'); hold on;
subplot(1,2,1);
plot(fq_small,abs(frf9));
xlabel('Frequency [Hz]');
title('FRF')
xlim([-500 500]);
subplot(1,2,2);
plot(fq,abs(fftshift(fft(Sy9))));
xlim([-500,500]);
xlabel('Frequency [Hz]');
title('Sy9')

figure;
subplot(5,2,1)
plot(fq_small,abs(FRF(Su,Sy,5)));
title('FRF period 1')
for i = 1:9
    subplot(5,2,i+1)
    u = sprintf('Su%d',i);
    y = sprintf('Sy%d',i);
    plot(fq_small,abs(FRF(data.(u),data.(y),5)));
    title(join(['FRF period ',num2str(i+1)]));
end

%%

% figure; hold on;
% subplot(1,2,1);
% plot(Su9);
% subplot(1,2,2);
% plot(fq,abs(fftshift(fft(Su9))));

load('Group09_Output2.Mat')


frf9 = FRF(Su9,Sy9,5);


figure('Name','output2'); hold on;
subplot(1,2,1);
plot(fq_small,abs(frf9));
title('FRF')
xlabel('Frequency [Hz]');
%xlim([-500 500]);
subplot(1,2,2);
plot(fq,abs(fftshift(fft(Sy9))));
xlim([-500,500]);
xlabel('Frequency [Hz]');
title('Sy9')

data = load('Group09_Output2.Mat');

figure;
subplot(10,2,1)
plot(Su);
title('Su')
subplot(10,2,2)
plot(Sy);
title('Sy')
for i = 1:9
    subplot(10,2,i*2+1)
    u = sprintf('Su%d',i);
    plot(data.(u));
    title(join(['Su',num2str(i)]));
    subplot(10,2,i*2+2)
    y = sprintf('Sy%d',i);
    plot(data.(y));
    title(join(['Sy',num2str(i)]));
end

figure;
subplot(10,2,1)
plot(fq,abs(fftshift(fft(Su))));
xlim([-510 510]);
title('Su')
subplot(10,2,2)
plot(fq,abs(fftshift(fft(Sy))));
xlim([-510 510]);
title('Sy')
for i = 1:9
    subplot(10,2,i*2+1)
    u = sprintf('Su%d',i);
    plot(fq,abs(fftshift(fft(data.(u)))));
    xlim([-510 510]);
    title(join(['Su',num2str(i)]));
    subplot(10,2,i*2+2)
    y = sprintf('Sy%d',i);
    plot(fq,abs(fftshift(fft(data.(y)))));
    xlim([-510 510]);
    title(join(['Sy',num2str(i)]));
end

figure;
subplot(5,2,1)
plot(fq_small,abs(FRF(Su,Sy,5)));
title('FRF period 1')
for i = 1:9
    subplot(5,2,i+1)
    u = sprintf('Su%d',i);
    y = sprintf('Sy%d',i);
    plot(fq_small,abs(FRF(data.(u),data.(y),5)));
    title(join(['FRF period ',num2str(i+1)]));
end

% figure; hold on;
% subplot(1,2,1);
% plot(Su9);
% subplot(1,2,2);
% plot(fq,abs(fftshift(fft(Su9))));


%%
load('Group09_Output3.Mat')


frf9 = FRF(Su9,Sy9,5);
% freqz(frf9);

figure('Name','output3'); hold on;
subplot(1,2,1);
plot(fq_small,abs(frf9));
title('FRF')
xlabel('Frequency [Hz]');
subplot(1,2,2);
plot(fq,abs(fftshift(fft(Sy9))));
xlim([-500 500]);
xlabel('Frequency [Hz]');
title('Sy9')

data = load('Group09_Output3.Mat');

figure;
subplot(10,2,1)
plot(Su);
title('Su')
subplot(10,2,2)
plot(Sy);
title('Sy')
for i = 1:9
    subplot(10,2,i*2+1)
    u = sprintf('Su%d',i);
    plot(data.(u));
    title(join(['Su',num2str(i)]));
    subplot(10,2,i*2+2)
    y = sprintf('Sy%d',i);
    plot(data.(y));
    title(join(['Sy',num2str(i)]));
end

figure;
subplot(10,2,1)
plot(fq,abs(fftshift(fft(Su))));
xlim([-510 510]);
title('Su')
subplot(10,2,2)
plot(fq,abs(fftshift(fft(Sy))));
xlim([-510 510]);
title('Sy')
for i = 1:9
    subplot(10,2,i*2+1)
    u = sprintf('Su%d',i);
    plot(fq,abs(fftshift(fft(data.(u)))));
    xlim([-510 510]);
    title(join(['Su',num2str(i)]));
    subplot(10,2,i*2+2)
    y = sprintf('Sy%d',i);
    plot(fq,abs(fftshift(fft(data.(y)))));
    xlim([-510 510]);
    title(join(['Sy',num2str(i)]));
end

figure;
subplot(5,2,1)
plot(fq_small,abs(FRF(Su,Sy,5)));
title('FRF period 1')
for i = 1:9
    subplot(5,2,i+1)
    u = sprintf('Su%d',i);
    y = sprintf('Sy%d',i);
    plot(fq_small,abs(FRF(data.(u),data.(y),5)));
    title(join(['FRF period ',num2str(i+1)]));
end

% figure; hold on;
% subplot(1,2,1);
% plot(Su9);
% subplot(1,2,2);
% plot(fq,abs(fftshift(fft(Su9))));


%%
load('Group09_Output4.mat')

frf1 = FRF(Su1,Sy1,0);
frf9 = FRF(Su9,Sy9,0);

[Hann] = hanning(length(Sy1),'periodic');

fq_S = -4000:4000-1;

figure;
subplot(1,3,1);
plot(fq_S',abs(fftshift(fft(Su1))));
title('Aperiodic');
subplot(1,3,2);
plot(fq_S,abs(fftshift(fft(Su9))));
title('Periodic');
subplot(1,3,3);
plot(fq_S,abs(fftshift(fft(Su9))) - abs(fftshift(fft(Su1))));
title('difference');

figure('Name','Noise outputs');
subplot(1,3,1);
plot(fq_S,abs(fftshift(fft(Sy1))));
title('Aperiodic');
xlim([-500 500]);
xlabel('Frequency [Hz]');
subplot(1,3,2);
plot(fq_S,abs(fftshift(fft(Sy9))));
title('Periodic');
xlim([-500 500]);
xlabel('Frequency [Hz]');
subplot(1,3,3);
plot(fq_S,abs(fftshift(fft(Sy1.*Hann))));
% plot(fq_S,abs(fftshift(fft(Sy9))) - abs(fftshift(fft(Sy1))));
title('difference');
title('Aperiodic times Hanning');
xlim([-500 500]);
xlabel('Frequency [Hz]');

figure; hold on;
subplot(1,3,1);
plot(fq,abs(frf9));
title('FRF periodic noise')
%xlim([-500 500]);
subplot(1,3,2);
plot(fq,abs(frf1));
title('FRF aperiodic noise')
subplot(1,3,3);

frf1 = FRF(Su1.*Hann,Sy1.*Hann,0);
plot(fq,abs(frf1));
title('FRF aperiodic noise times Hann')

figure; hold on;
plot(fq,abs(frf9));
xlabel('Frequency [Hz]');
xlim([-500 500]);
title('FRF periodic noise')
%xlim([-500 500]);
figure; hold on;
plot(fq,abs(frf1));
xlabel('Frequency [Hz]');
xlim([-500 500]);
title('FRF aperiodic noise')
figure; hold on;
[Hann] = hanning(length(Sy1),'periodic');
frf1 = FRF(Su1.*Hann,Sy1.*Hann,0);
plot(fq,abs(frf1));
xlabel('Frequency [Hz]');
xlim([-500 500]);
title('FRF aperiodic noise times Hann')

% y = randn(1,N);
% y = 0.1*y/rms(y);
% 
% save('Group09_Input4.Mat','y')
