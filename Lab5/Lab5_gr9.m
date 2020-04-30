clear;close all;clc;

%%
% 
% Exact model
% 
% 

% Continuous time Butterworth filter
% Filter order : 2
% Cutoff angular frequency : [0.01,0.2]

[B_tilde,A_tilde] = butter(1,[0.01,0.2],'s');
B_tilde = B_tilde/A_tilde(end);
A_tilde = A_tilde/A_tilde(end);

w_start = 0.001;
w_stop = 2;

N = 500;            % Number of frequencies

W = linspace(w_start,w_stop,N)';

% Exact frequency response of the reference filter
G_0 = freqs(B_tilde,A_tilde,W);
freqs(B_tilde,A_tilde,W)

figure
semilogx(db(G_0),'+')
%% Circular zero mean white noise
sigma = 0.001*0;          % Standard deviation
noise = sigma*(randn(N,1) + 1i*randn(N,1))/sqrt(2);

% Measured frequncy response
G_m = G_0 + noise;

% Rewrite the transfer funciton
nB = length(B_tilde)-1;
nA = length(A_tilde)-1;
s = 1i*W;
%%
maxpow = max(nB,nA);
s_all = repmat(s,[1,maxpow+1]).^(repmat((0:maxpow),[N,1]));

%%
% $B(s,a) = b_{2}s^{2}+b_{1}s+b_{0}$
%B = (B_tilde).*s_all/A_tilde(1);
%B = sum(B,1);
B0 = polyval(B_tilde,s);

%%
% $A'(s,a) = a_{2}s^{2}+a_{1}s$
%A_prime = A_tilde(2:end)/A_tilde(1);
%A_prime = A_prime.'.*s(2:end,:);
%A_prime = sum(A_prime,1);
A0 = polyval(A_tilde,s);    % a_2*s^2 + a_1*s + a_0

% True parameters
theta_true = [B_tilde A_tilde(1:end-1)]; %b2 b1 b0 a2 a1
%theta_true = fliplr(theta_true).'; 


%% Plots
figure;
Wlog = logspace(-3,1,N);
loglog(Wlog,abs(G_0));
hold on;
loglog(Wlog,abs(G_m));
grid on;
legend('Exact frequency response','Measured frequncy response')
xlabel('Frequency (rad/s)');
ylabel('Magnitude');

figure;
plot(W,db(G_0));
hold on;
plot(W,db(G_m));
grid on;
legend('Exact frequency response','Measured frequncy response');


%%
% 
%  Implemenation of the Levy estimators
%  
% 

% Cost function vectorization
%s = s.'; 
% e_Levy = G_m + G_m.*(A_prime.') - B.';    % For perfect reconstruction
e_Levy = G_m ;

% Jacobian . Theta: %b2 b1 b0 a2 a1
J_Levy = [-s_all(:,3) -s_all(:,2) -s_all(:,1) G_m.*s_all(:,3) G_m.*s_all(:,2)];
%%
% Real/imaginary part separation
e_Levy_IR = [real(e_Levy) ; imag(e_Levy)]; 
J_Levy_IR = [real(J_Levy) ; imag(J_Levy)];

% Levy theta estimation
theta_Levy = -J_Levy_IR\e_Levy_IR;

%%
Thetaest_B = theta_Levy(1:nB+1);
Thetaest_Ap = theta_Levy(nB+2:end);

GestLevy = freqs(Thetaest_B,[Thetaest_Ap; 1],W);
%%
figure
plot(W,db(GestLevy),'o',W,db(G_m),'*',W,db(G_0),'+')

keyboard
%%
disp('Parameters comparison');
disp(join(['True theta:           ',num2str(theta_true.')]));
disp(join(['Levy estimated theta: ',num2str(theta_Levy(1:end-1).')]));
% disp(join(['Difference:           ',num2str(theta_true.' - theta_Levy(1:end-1).')]));

%%
% 
%  Implemenation of the Sanathanan estimator
%  
% 

% Iteration index
l_max = 8;

% Computing B with new parameters
B = flip(theta_Levy(1:3).').*s(:,1:3);
B = sum(B,2);

% Computing A with new parameters
A_prime = theta_Levy(4:5);
A_prime = flip(A_prime.').*s(:,2:3);
A_prime = sum(A_prime,2);

% Cost function vectorization
e_San = (G_m + G_m.*(A_prime) - B)./vecnorm(1+A_prime,2,2);

% Jacobian
J_San = J_Levy./vecnorm(1+A_prime.',2,2);

% Real/imaginary part separation
e_San_IR = [real(e_San) ; imag(e_San)]; 
J_San_IR = [real(J_San);imag(J_San)];

theta_San = theta_Levy; 

% Sanathanan estimation
for l = 1:l_max
    
    % Parameters computation
    theta_San = J_San_IR\e_San_IR;
    theta_San = theta_San(1:end-1);
    
    % Computing B with new parameters
    B = flip(theta_San(1:3).').*s(:,1:3);
    B = sum(B,2);
    
    % Computing A with new parameters
    A_prime_new = theta_San(4:5);
    A_prime_new = flip(A_prime_new.').*s(:,2:3);
    A_prime_new = sum(A_prime_new,2);
    
    % Updating the cost function and the Jacobian
    e_San = (G_m + G_m.*(A_prime_new) - B)./vecnorm(1+A_prime.',2,2);
    J_San = J_Levy./vecnorm(1+A_prime.',2,2);
    
    e_San_IR = [real(e_San) ; imag(e_San)]; 
    J_San_IR = [real(J_San);imag(J_San)];
    
    % Sacing the previous A_prime
    A_prime = A_prime_new.';
end

% Parameters computation
theta_San = J_San_IR\e_San_IR;
theta_San = theta_San(1:end-1);

disp(join(['Sanathanan estimated theta: ',num2str(theta_San.')]));
%%
% 
%  Implemenating Gauss-Newton based least squares estimation
%  
% 


