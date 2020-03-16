P=tf(1,[1 1])
bode(P)
% design PI controller for crossover =10, PM=50
a=10/tan(44*pi/180)
C=tf([1 a],[1 0])
bode(P*C)
k=1/bode(P*C,10)
C=tf([1 a],[1 0])*k
%evaluate feedback system
bode(P*C) % loop tf
bode(feedback(P*C,1),P) % closed-loop sys freq response, check BW
step(feedback(P*C,1),P) % response to steps in reference
bode(feedback(P,C),P)   % freq response to disturbance
step(feedback(P,C),P)   % response to step disturbance

% design PI controller for crossover =10, PM=70
a=10/tan(64*pi/180)
C2=tf([1 a],[1 0])
k=1/bode(P*C2,10)
C2=tf([1 a],[1 0])*k
%evaluate feedback system
bode(P*C,P*C2,P) % loop tf
bode(feedback(P*C,1),feedback(P*C2,1),P) % closed-loop sys freq response, check BW
step(feedback(P*C,1),feedback(P*C2,1),P) % closed-loop sys step response
bode(feedback(P,C),feedback(P,C2),P) % freq response to disturbance
step(feedback(P,C),feedback(P,C2),P)  % response to step disturbance
uc=tf([1 -1],[1 1])
nyquist(P*C,P*C2,P,uc) % nyquist plots
step(feedback(C,P),feedback(C2,P))  % control input for steps in reference