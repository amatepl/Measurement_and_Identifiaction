function H = HFunction_time(u,y)

u = sum(u,2)/size(u,2);
y = sum(y,2)/size(u,2);

H = fft(y)./fft(u);
H = fftshift(H);