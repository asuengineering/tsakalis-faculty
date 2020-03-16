% SVD application in the modeling of a noisy oscillatory signal 
%     as the output of an autoregressive model:
%     y(n)=[y(n-1),y(n-2),...]*q

% Define the signals
t=(0:.01:10)';          % time
n=rand(size(t))-.5;     % noise
y=sin(10*t+1)+.02*n;    % measurement
N0=10:length(t);        % fitting window

% Optional Filtering for frequency-weighted fit
%F=tf(.01,[1 -.99],.01);y=lsim(F,y); 

% Form the regressor by taking lags of the output
W=[y(N0-1),y(N0-2),y(N0-3),y(N0-4),y(N0-5),y(N0-6),y(N0-7),y(N0-8),y(N0-9)];

% least squares fit
q=W\(y(N0))
plot(t,y,t(N0),W*q,t(N0),W*q-y(N0)); pause    % check the fit

% Autoregressive transfer function: resonance at the oscillation frequency
g=tf(1,[1 -q'],.01)
bode(g)                 % check the t.f.

% Model order: How many columns of W do you need?
s=svd(W)            % svd of regressor
s2=svd(W'*W)         % svd of gramian

% Questions:
% 1. What is the relationship between s and s2? How many lags do you need
% in the model?   
% 2. The singular values of W appear to reach a floor related to the noise. 
%    Derive this value analytically and verify with an example.
% 3. What is ithe effect of the noise amplitude?
% 4. What happens when the signal is composed of two frequencies, say 10
% and 2?
