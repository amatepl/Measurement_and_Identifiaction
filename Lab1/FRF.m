%%
%
% Compute the Frequency Response Function
%
%

function H = FRF(Su,Sy)

% FFT
U = fftshift(fft(Su));
Y = fftshift(fft(Sy));

%%
%
% In order to avoid the division by 0 and by that making shure we compute
% the FRF only for the frequency band of interest the computation is done
% using the input values higher (in the frequency domain) than a given 
% threshold and the corresponding output.
%

% FRF computation
threshold = 6; 

H = Y(U >threshold)./U(U > threshold);