P=tf(1,[1 1])
w=10
bode(P)
C=tf(1,[1 0])
bode(P*C)
[m,p]=bode(P*C,w)
p+120
tz=tan(54*pi/180)/w
1/tz
C=tf([tz 1],[1 0])
bode(P*C)
grid
[m,p]=bode(P*C,w)
C=tf([tz 1],[1 0])/m
bode(P*C)
bode(feedback(P*C,1))
step(feedback(P*C,1))
step(feedback(P,C))
C2=tf(10*[1 1],[1 0])
bode(feedback(P*C,1),feedback(P*C2,1))
bode(feedback(1,P*C),feedback(1,P*C2))
step(feedback(P*C,1),feedback(P*C2,1))
step(feedback(P,C),feedback(P,C2))
bode((P*C),(P*C2))
T=.03;step(feedback(P*C,1),feedback(c2d(P,T,'zoh')*c2d(C,T,'tustin'),1))
T=.03;bode((P*C),(c2d(P,T,'zoh')*c2d(C,T,'tustin')))
T=.03;step(feedback(P*C,1),feedback(c2d(P,T,'zoh')*c2d(C,T,'zoh'),1))
T=.06;step(feedback(P*C,1),feedback(c2d(P,T,'zoh')*c2d(C,T,'tustin'),1))
T=.6;step(feedback(P*C,1),feedback(c2d(P,T,'zoh')*c2d(C,T,'tustin'),1))
T=.003;step(feedback(P*C,1),feedback(c2d(P,T,'zoh')*c2d(C,T,'tustin'),1))
T=.003;c2d(C,T,'tustin')
T=.03;c2d(C,T,'tustin')
T=.000003;c2d(C,T,'tustin')
tz=tan(63*pi/180)/w
C2=tf([tz 1],[1 0])
[m,p]=bode(P*C2,w)
C2=tf([tz 1],[1 0])/m
T=.03;step(feedback(P*C,1),feedback(c2d(P,T,'zoh')*c2d(C2,T,'tustin'),1))
bode((P*C),(P*C2))
T=.03;bode((P*C),(c2d(P,T,'zoh')*c2d(C2,T,'tustin')))
