P=tf(1,[1 0])
wc=1,PM=60,
bode(P)
[m,p]=bode(P,wc);
z_angle=-((p-90)-(-180)-PM)
T=max(0,tan(z_angle*pi/180)/wc)
C=tf([T 1],[1 0])
k=1./bode(P*C,wc)
C=tf([T 1],[1 0])*k
step(feedback(P*C,1))
bode(feedback(P*C,1))
pi*60/90
ts=.2;step(feedback(P*C,1),feedback(c2d(P,ts)*c2d(C,ts,'tustin'),1))
ts=1;step(feedback(P*C,1),feedback(c2d(P,ts)*c2d(C,ts,'tustin'),1))
ts=2;step(feedback(P*C,1),feedback(c2d(P,ts)*c2d(C,ts,'tustin'),1))
ts=2.1;step(feedback(P*C,1),feedback(c2d(P,ts)*c2d(C,ts,'tustin'),1))
ts=2.2;step(feedback(P*C,1),feedback(c2d(P,ts)*c2d(C,ts,'tustin'),1))
ts=2.3;step(feedback(P*C,1),feedback(c2d(P,ts)*c2d(C,ts,'tustin'),1))
ts=2.5;step(feedback(P*C,1),feedback(c2d(P,ts)*c2d(C,ts,'tustin'),1))

