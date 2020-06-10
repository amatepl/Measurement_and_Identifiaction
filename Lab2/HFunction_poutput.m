function H = HFunction_poutput(u,y)

U = fft(u);
Y = fft(y);

S_YY = sum(Y.*(conj(Y)),2)/size(U,2);
S_UY = sum(U.*(conj(Y)),2)/size(U,2);

H = S_YY./S_UY;
H = fftshift(H);