function H = HFunction_freq(u,y)

U = sum(fft(u,[],1),2)/size(u,2);
Y = sum(fft(y,[],1),2)/size(u,2);

H = Y./U;
H = fftshift(H);
