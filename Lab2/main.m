clear; close all;clc;
%% Multisine generation
% Parameters
N1 = 4000;
V_RMS = 0.5;
K_min = 4/4;
K_max = 1000/4;
K = K_min:K_max;

% Random phase multisine
phi = 2*pi*rand([K_min,K_max]);

% Frequency domain construction
X = zeros(N1,1);
X(K_min+1:K_max+1) = exp(1j*phi);
multisine_rand = N1*real(ifft(X));
multisine_rand = multisine_rand*V_RMS/rms(multisine_rand);

noise = randn([1,160000]);
noise = noise*V_RMS/rms(noise);


% Schroeder phase multisine
phi = K.*(K+1)*pi/length(K);

% Frequency domain construction
X = zeros(N1,1);
X(K_min+1:K_max+1) = exp(1j*phi);
multisine_scho = N1*real(ifft(X));
multisine_scho = multisine_scho*V_RMS/rms(multisine_scho);
fq = -N1/2:N1/2 - 1;
% fq = -N1*2:4:N*2 - 1;
fq = fq*4;
figure; hold on;
subplot(1,2,1);
plot(fq,abs(fftshift(fft(multisine_scho))));
subplot(1,2,2);
plot(multisine_scho);

%%

[u , y] =  ReadDataLab2(N1,40,32,'out_Group9_InputPeriodick.mat');

[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_time);

figure('Name','Time');
subplot(1,2,1);
plot(db(H(:,1)),'b.');
hold on;
plot(db(H(:,2)),'r.');

subplot(1,2,2);
plot((stdH(:,1)),'bo');
hold on;
plot((stdH(:,4)),'ro');

[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_freq);

figure('Name','Freq');
subplot(1,2,1);
plot(db(H(:,1)),'b.');
hold on;
plot(db(H(:,5)),'r.');

subplot(1,2,2);
plot((stdH(:,1)),'bo');
hold on;
plot((stdH(:,4)),'ro');

[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_FRF);

figure('Name','FRF');
subplot(1,2,1);
plot(db(H(:,1)),'b.');
hold on;
plot(db(H(:,5)),'r.');

subplot(1,2,2);
plot((stdH(:,1)),'bo');
hold on;
plot((stdH(:,4)),'ro');

[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_pinput);

figure('Name','Power 1');
subplot(1,2,1);
plot(db(H(:,1)),'b.');
hold on;
plot(db(H(:,5)),'r.');

subplot(1,2,2);
plot((stdH(:,1)),'bo');
hold on;
plot((stdH(:,4)),'ro');

[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_poutput);

figure('Name','Power 2');
subplot(1,2,1);
plot(db(H(:,1)),'b.');
hold on;
plot(db(H(:,5)),'r.');

subplot(1,2,2);
plot((stdH(:,1)),'bo');
hold on;
plot((stdH(:,4)),'ro');

%% Plots Aperiodic 

% [u , y] =  ReadDataLab2(N1*40,40,32,'out_Group9_InputAperiodic.mat');

data = load('out_Group9_InputAperiodic.mat');
u = data.Su;
y = data.Sy;

u = reshape(u,[],40);
y = reshape(y,[],40); 

u = u(:,8:end);
y = y(:,8:end);

[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_time);

figure('Name','Time');
subplot(1,2,1);
plot(db(H(:,1)),'b.');
hold on;
plot(db(H(:,5)),'r.');

subplot(1,2,2);
plot((stdH(:,1)),'bo');
hold on;
plot((stdH(:,4)),'ro');

[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_freq);

figure('Name','Freq');
subplot(1,2,1);
plot(db(H(:,1)),'b.');
hold on;
plot(db(H(:,5)),'r.');

subplot(1,2,2);
plot((stdH(:,1)),'bo');
hold on;
plot((stdH(:,4)),'ro');

[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_FRF);

figure('Name','FRF');
subplot(1,2,1);
plot(db(H(:,1)),'b.');
hold on;
plot(db(H(:,5)),'r.');

subplot(1,2,2);
plot((stdH(:,1)),'bo');
hold on;
plot((stdH(:,4)),'ro');

[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_pinput);

figure('Name','Power 1');
subplot(1,2,1);
plot(db(H(:,1)),'b.');
hold on;
plot(db(H(:,5)),'r.');

subplot(1,2,2);
plot((stdH(:,1)),'bo');
hold on;
plot((stdH(:,4)),'ro');

[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_poutput);

figure('Name','Power 2');
subplot(1,2,1);
plot(db(H(:,1)),'b.');
hold on;
plot(db(H(:,5)),'r.');

subplot(1,2,2);
plot((stdH(:,1)),'bo');
hold on;
plot((stdH(:,4)),'ro');
