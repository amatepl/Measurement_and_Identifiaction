close all;clear;clc;

% Sampling frequency
fs = 1e4;
% Exited frequency
F = 100;
% Number of samples
N = 4096;
% time
t = (0:4096-1)/fs;
% Frequency axis
f=(0:N-1)*fs/N;


first_signal = sin(2*pi*F*t');

figure;
plot(t',first_signal)
xlabel('Time (s)');
ylabel('Amplitude');

figure;
plot(f,db((fft(first_signal))),'o')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

%%
% Set the number of periods for leakage correction
P = 8;
F = 5*P*fs/N;

first_signal = sin(2*pi*F*t');

figure;
plot(f,db((fft(first_signal))),'o')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

%%
% Second signal
A = logspace(log10(0.1),log10(1.1),10);

second_signal = first_signal.*ones(length(first_signal),10);
second_signal = second_signal.*A;
fft_second_signal = fft(second_signal,[],1);
second_signal = reshape(second_signal,[],1);

figure;
plot(second_signal);
xlabel('Time (s)');
ylabel('Amplitude');

figure;
plot(f,db(fft_second_signal),'o');
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

%%
load('measGroup9_Input1.mat');

f = f(1:4096/P);

figure('Name','First signal');
plot(f,db(fft(Su(end-4096/P +1:end))),'o',f,db(fft(Sy(end-4096/P +1:end))),'o');
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Su','Sy')
%%
load('measGroup9_Input2.mat');

figure('Name','First signal');
plot(f,db(fft(Su(end-4096/P +1:end))),'o',f,db(fft(Sy(end-4096/P +1:end))),'o');
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Su','Sy')

figure('Name','First signal');
plot(f,db(fft(Su(end-4096/P +1:end))),'o');
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Su')

figure('Name','First signal');
plot(f,db(fft(Sy(end-4096/P +1:end))),'o');
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Sy')


%%
load('measGroup9_Input3.mat');

% figure('Name','First signal');
% for i = 1:10
%     Su = reshape(Su,4096,10);
%     Sy = reshape(Sy,4096,10);
%     subplot(5,2,i);
%     plot(f,db(fft(Su(end-4096/P +1:end,i))),'o',f,db(fft(Sy(end-4096/P +1:end,i))),'o');
%     title(num2str(i))
%     xlabel('Frequency (Hz)')
%     ylabel('Magnitude (dB)')
%     legend('Su','Sy')
% end

Pin = zeros(10,1);
Pout = zeros(10,1);

for i = 1:10
    figure('Name','First signal');
    Su = reshape(Su,4096,10);
    Sy = reshape(Sy,4096,10);
    plot(f,db(fft(Su(end-4096/P +1:end,i))),'o',f,db(fft(Sy(end-4096/P +1:end,i))),'o');
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    legend('Su','Sy')
    
    Pin(i) = rms(Su(end-4096/P +1:end,i))^2;
    Pout(i) = rms(Sy(end-4096/P +1:end,i))^2;
    
%     figure('Name','First signal');
%     Su = reshape(Su,4096,10);
%     Sy = reshape(Sy,4096,10);
%     
%     loglog(db((Su(end-4096/P +1:end,i))),db((Sy(end-4096/P +1:end,i))));
%     xlabel('Su (dB)')
%     ylabel('Sy (dB)')
end

Pidbm = 10*log10(Pin)+30;
Podbm = 10*log10(Pout)+30;

K = (Pout(2) - Pout(1))/(Pin(2)-Pin(1));

Kdbm = (Podbm(2) - Podbm(1))/(Pidbm(2)-Pidbm(1));

figure;
datacursormode on;
% plot(Pin - Pin(1),Pout - Pout(1))
hold on;
% plot(Pin,K*Pin+Pout(1)-Pin(1))

plot(10*log10(Pin)+30,10*log10(Pout)+30);
hold on;
plot(Pidbm,Kdbm*Pidbm + Podbm(1) - Kdbm*Pidbm(1))
xline(16.244250726673570)
legend('Studied system','Linear system');
xlabel('Input power (dBm)');
ylabel('Output power (dBm)');
% loglog(10*log10(Pin)+30,(10*log10(K*Pin+Pout(1)-Pin(1))+30));
% plot(10*log10(Pin)+30,10*log10(Pin)+30);
% set(gca, 'XScale', 'log', 'YScale', 'log');
