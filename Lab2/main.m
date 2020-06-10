clear; close all;clc;
%% Multisine generation
% Parameters
N1 = 4000;          % Number of samples
V_RMS = 0.5;        % RMS value of the input signal
K_min = 4/4;        % Minimum frequency bin. (The sampling frequency is 4 times larger than N)
K_max = 1000/4;     % Maximum frequency bin
K = K_min:K_max;

% Random phase multisine
phi = 2*pi*rand([K_min,K_max]);     % Random phase

% Frequency domain construction
X = zeros(N1,1);
X(K_min+1:K_max+1) = exp(1j*phi);   % Signal built in the frequency domain
multisine_rand = N1*real(ifft(X));  % Time domain
multisine_rand = multisine_rand*V_RMS/rms(multisine_rand);  % Setting the RMS value

% Noise
noise = randn([1,160000]);
noise = noise*V_RMS/rms(noise);     % Setting the RMS value

% Schroeder phase multisine
phi = K.*(K+1)*pi/length(K);        % Schroeder phase

% Frequency domain construction
X = zeros(N1,1);
X(K_min+1:K_max+1) = exp(1j*phi);   % Signal built in the frequency domain
multisine_scho = N1*real(ifft(X));  % Time domain
multisine_scho = multisine_scho*V_RMS/rms(multisine_scho);   % Setting the RMS value

fq = -N1/2:N1/2 - 1;
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

%%
[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_time);
figure('Name','Time');
plot(fq,db(H(:,1)),'b.');
hold on;
plot(fq,db(H(:,5)),'r.');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
legend('2 repetitions','32 repetitions');
[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_freq);
figure('Name','Time');
plot(fq,db(H(:,1)),'b.');
hold on;
plot(fq,db(H(:,5)),'r.');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
legend('2 repetitions','32 repetitions');
[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_FRF);
figure('Name','Time');
plot(fq,db(H(:,1)),'b.');
hold on;
plot(fq,db(H(:,5)),'r.');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
legend('2 repetitions','32 repetitions');
[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_pinput);
figure('Name','Time');
plot(fq,db(H(:,1)),'b.');
hold on;
plot(fq,db(H(:,5)),'r.');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
legend('2 repetitions','32 repetitions');
[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_poutput);
figure('Name','Time');
plot(fq,db(H(:,1)),'b.');
hold on;
plot(fq,db(H(:,5)),'r.');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
legend('2 repetitions','32 repetitions');

%%
[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_time);
figure('Name','Time');
plot(fq,db(stdH(:,1)),'b.');
hold on;
plot(fq,db(stdH(:,4)),'r.');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
legend('2 repetitions','32 repetitions');
[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_freq);
figure('Name','Time');
plot(fq,db(stdH(:,1)),'b.');
hold on;
plot(fq,db(stdH(:,4)),'r.');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
legend('2 repetitions','32 repetitions');
[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_FRF);
figure('Name','Time');
plot(fq,db(stdH(:,1)),'b.');
hold on;
plot(fq,db(stdH(:,4)),'r.');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
legend('2 repetitions','32 repetitions');
[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_pinput);
figure('Name','Time');
plot(fq,db(stdH(:,1)),'b.');
hold on;
plot(fq,db(stdH(:,4)),'r.');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
legend('2 repetitions','32 repetitions');
[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_poutput);
figure('Name','Time');
plot(fq,db(stdH(:,1)),'b.');
hold on;
plot(fq,db(stdH(:,4)),'r.');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
legend('2 repetitions','32 repetitions');
%%
[H, stdH] = TransferFunc(u, y, [2,4,8,16,32], @HFunction_time);
figure('Name','Time');
for i = 1:5
    subplot(2,3,i)
    plot(fq,db(H(:,i)),'.');
    hold on;
    ylabel('Magnitude (dB)');
    xlabel('Frequency (Hz)');
    title(join([num2str(2^(i)),' repetitions']));
    %legend(join([num2str(2^(i)),' repetitions']));
end


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


