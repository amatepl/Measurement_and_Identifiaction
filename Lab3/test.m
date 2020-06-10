
N = 1000;
fs = 1000;
t = (1:N)/fs;
x = sin(2*pi*4*t) + sin(2*pi*11*t);
y = x - (x.^3)/2 - (x.^4)/4;
dft = fftshift(fft(y));

f = -500:499;
figure;
% plot(y)
stem(f,abs(dft),'o');