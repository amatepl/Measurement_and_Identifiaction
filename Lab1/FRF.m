%%
%
% Compute the Frequency Response Function
%
%

function H = FRF(Su,Sy,threshold)

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
% threshold = 0.00; 

H = Y(abs(U) >threshold)./U(abs(U) > threshold);

% threshold = 5.0; 
% H(abs(U) <threshold)=0;