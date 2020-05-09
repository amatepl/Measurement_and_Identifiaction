function H = HFunction_time(u,y)

u = sum(u,2)/size(u,2);
y = sum(y,2)/size(u,2);

u = mean(u,2);
y = mean(y,2);

figure;
plot(u);
hold on;
plot(y);

H = fft(y)./fft(u);
H = fftshift(H);
% H = H/size(u,2);