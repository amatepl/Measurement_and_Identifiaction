clear, close all, clc;

% DFT of 3 periods of a cosine.

nPeriods = 3;   % Number of periods in N samples
N = 1000;       % Number of samples

omega = 2*pi*nPeriods/N;
phi = randi(N,[1, 3]);  % phase

t = 1:N;
fq = -N/2:N/2-1;
cosinus = cos(omega*t + phi.');
% cosinus = cos(omega*t);
cosinus = sum(cosinus,1);

cos_fft = fftshift(fft(cosinus,[],2));
tresh = 10^(-10);
cos_fft(abs(cos_fft)<tresh) =0;

%% Plots
figure; hold on;
plot(t,cos(omega*t + phi.'));
% plot(t,cos(omega*t));
plot(cosinus);
legend(join(['\phi = ',num2str(phi(1))]),join(['\phi = ',num2str(phi(2))]),join(['\phi = ',num2str(phi(3))]),'Sum of the sequence');
title('Cosine sequence with a randomly selected phase');

figure; hold on;
plot(fq,abs(cos_fft));
title('DFT of the sum of cosinus sequence')

figure; hold on;
subplot(2,1,1);
stem(fq,abs(cos_fft));
title('Frequency axis in bins');
subplot(2,1,2);
stem(fq,angle(cos_fft)/pi);
xlabel('Frequency [bins]');

fs = 100;
figure; hold on;
subplot(2,1,1);
stem(fq*fs/N,abs(cos_fft));
title('Frequency axis in Hz');
subplot(2,1,2);
stem(fq*fs/N,angle(cos_fft));
xlabel('Frequency [Hz]');



% %%
% fs = 100;
% t = 0:1/fs:1-1/fs;
% x = cos(2*pi*15*t - pi/4) - sin(2*pi*40*t);
% 
% y = fft(x);
% z = fftshift(y);
% 
% ly = length(y);
% f = (-ly/2:ly/2-1)/ly*fs;
% 
% stem(f,abs(z))
% xlabel 'Frequency (Hz)'
% ylabel '|y|'
% grid
% 
% tol = 1e-6;
% z(abs(z) < tol) = 0;
% 
% theta = angle(z);
% 
% stem(f,theta/pi)
% xlabel 'Frequency (Hz)'
% ylabel 'Phase / \pi'
% grid