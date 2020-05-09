function H = HFunction_pinput(u,y)

U = fft(u);
Y = fft(y);

S_YU = sum(Y.*(conj(U)),2)/size(U,2);
S_UU = sum(U.*(conj(U)),2)/size(U,2);

H = S_YU./S_UU;