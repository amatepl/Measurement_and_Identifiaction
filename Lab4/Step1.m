close all; clear;clc;

% 4.3.1 Build the discrite time test

filterOrder = 2;
ripple = 3;             % 3dB
range = [0.3 0.6];

[num,denom] = cheby2(filterOrder,ripple,range);

F = tf(num,denom,1);

[F,t] = impulse(F);

G = F(1:25);


% 4.3.2 Generate the excitation and the output signals.

Ne = size(G,1)*10;
ue = randn([Ne,1]);
uv = randn([Ne,1]);

ye = filter(num,denom,ue);
yv = filter(num,denom,uv);

% Noise

SNR = 6;
vare = var(ye)/(10^(SNR/10));
sde = sqrt(vare);
noisee = sde*randn([Ne,1]);

varv = var(yv)/(10^(SNR/10));
sdv = sqrt(varv);
noisev = sdv*randn([Ne,1]);

ye = ye + noisee;
yv = yv + noisev;

np = 10;
He = toeplitz(ue,ue(1:np));
Hv = toeplitz(uv,uv(1:np));

thetae = He\ye;
thetav = Hv\yv;