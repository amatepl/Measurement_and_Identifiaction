close all; clear;clc;



%%% System Parameters ----------------------------------------------------
filterOrder = 2;
ripple = 3;             % 3dB
range = [0.3 0.6];      % initial = [0.3 0.6]. Decrease for part 4.3.5
truncation = 38;
Ne = truncation*10;
SNR = 6;            % dB
np = 100;           % Maximum number of parameters
Nruns = 100;        % Maximum number of experiment runs
%%%-----------------------------------------------------------------------

% 4.3.1 Build the discrite time test

% Chebycheff filter creation

[num,denom] = cheby2(filterOrder,ripple,range);
figure;hold on;
freqz(num,denom);
F = tf(num,denom,1);
[F,t] = impulse(F);

G = F(1:truncation);                % G is the model

figure; hold on;
plot(F);
xlabel('n');
ylabel('F');

% freq(G);
%%

V = zeros([np,Nruns]);      % Least squares estimation cost function
VAIC = zeros([np,Nruns]);
Vval = zeros([np,Nruns]);   % Least squares validaton cost function

Vmin = zeros([Nruns,1]);      % Least squares estimation cost function
VAICmin = zeros([Nruns,1]);
Vvalmin = zeros([Nruns,1]);   % Least squares validaton cost function

for run = 1:Nruns
    
    % 4.3.2 Generate the excitation and the output signals.
    ue = randn([Ne,1]);     % Dataset for the estimation of the model
    uv = randn([Ne,1]);     % Dataset for the validation of the model

    ye = filter(num,denom,ue);  % Noiseless output of the system
    yv = filter(num,denom,uv);
    
    % Noise

    % Required SNR values for the tests are 6 and 26 db         
    vare = var(ye)/(10^(SNR/10));
    sde = sqrt(vare);
    noisee = sde*randn([Ne,1]);

    varv = var(yv)/(10^(SNR/10));
    sdv = sqrt(varv);
    noisev = sdv*randn([Ne,1]);

    ye = ye + noisee;       % Output of the system
    yv = yv + noisev;
    
    % 4.3.3 Implement the lest squares for varying n

    % Observation matrix
    He = toeplitz(ue,ue(1:np));
    Hv = toeplitz(uv,uv(1:np));

    
    
    for n = 1:np
        % Observation matrix
        Hel = He(:,1:n);
        Hvl = Hv(:,1:n);

        % Parameters
        theta = Hel\ye;     

        % Cost functions
        V(n,run) = (norm(ye - Hel*theta)^2)/(Ne*var(ue));
        VAIC(n,run) = V(n,run)*(1+2*n/Ne); 
        Vval(n,run) = (norm(yv - Hvl*theta)^2)/(Ne*var(uv));
    end

    % Minimums
    [~,Vmin(run)] = min(V(:,run));
    [~,VAICmin(run)] = min(VAIC(:,run));
    [~,Vvalmin(run)] = min(Vval(:,run));
    
end



% figure;hold on;
% plot(V);
% ylabel('V_{LS}');
% xlabel('parameter order n');
% title('Normalized least squares estimation cost function');

% 4.3.4 Select the optimal model order using AIC and Validation

% figure;hold on;
% plot(VAIC);
% ylabel('V_{AIC}');
% xlabel('parameter order n');
% title('AKAIKE information criteria');

% figure;hold on;
% plot(Vval);
% ylabel('V_{LS,val}');
% xlabel('parameter order n');
% title('Normalized least squares validation cost function');
tit = join(['Normalized least squares cost functions for SNR = ',num2str(SNR),' dB']);

figure;hold on;
plot(V(:,1));
plot(VAIC(:,1));
plot(Vval(:,1));
grid on;
legend('V_{LS}','V_{AIC}','V_{val}');
ylabel('Cost');
xlabel('parameter order n');
title(tit);

% 4.3.5 Explore the effects of SNR and bandwidth.
% Set different values of the SNR above

% 4.3.6 Determine the robustnes of the model order estimation
% figure;hold on;
% hist(Vmin,1:np);
% xlabel('Number of runs');
% % title(join(['Estimated model order for V_{LS} with ',num2str(Nruns),' runs']));
% 
% figure;hold on;
% hist(VAICmin,1:np);
% xlabel('Number of runs');
% % title(join(['Estimated model order for V_{AIC} with ',num2str(Nruns),' runs']));
% 
% figure;hold on;
% hist(Vvalmin,1:np);
% xlabel('Number of runs');
% % title(join(['Estimated model order for V_{val} with ',num2str(Nruns),' runs']));

figure('Position',[500 500 2000 500]);hold on;
set(gca,'DataAspectRatio',[1,1,1])
title(join(['SNR = ',num2str(SNR)]));
subplot(1,3,1);
hist(Vmin,1:np);
xlabel('Number of runs');
title("V_{LS}");
subplot(1,3,2);
hist(VAICmin,1:np);
xlabel('Number of runs');
title('V_{AIC}');
subplot(1,3,3);
hist(Vvalmin,1:np);
xlabel('Number of runs');
title('V_{val}');


