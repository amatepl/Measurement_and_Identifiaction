function H = HFunction_FRF(u,y)

H = fft(y,[],1)./fft(u,[],1);
H = sum(H,2)/size(H,2);