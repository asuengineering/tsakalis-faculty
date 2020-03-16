%-- 4/4/07  7:26 PM --%
P=tf(1,[1 1]);P.iodelay=0.1;
Pp=pade(P,2)
step(Pp)
L=0.09
R=(0.301-0.101)/(0.458-0.206)
Kp=1.2/R/L;Ki=0.6/R/L/L;Kd=0.6/R;
s=tf([1 0],1);
C=Kp+Ki/s+Kd*s/(.01*s+1)
step(feedback(C*Pp,1))
bode(C*Pp)
grid
bode(C*Pp,C)
L=0.2
Kp=1.2/R/L;Ki=0.6/R/L/L;Kd=0.6/R;C=Kp+Ki/s+Kd*s/(.01*s+1)
bode(C*Pp,C)
grid
step(feedback(C*Pp,1))
L=0.4
Kp=1.2/R/L;Ki=0.6/R/L/L;Kd=0.6/R;C=Kp+Ki/s+Kd*s/(.01*s+1)
bode(C*Pp)
step(feedback(C*Pp,1))
R=R*2
L=0.2
Kp=1.2/R/L;Ki=0.6/R/L/L;Kd=0.6/R;C=Kp+Ki/s+Kd*s/(.01*s+1)
step(feedback(C*Pp,1))
step(feedback(10*Pp,1))
step(feedback(15*Pp,1))
step(feedback(16*Pp,1))
step(feedback(18*Pp,1))
step(feedback(17*Pp,1))
step(feedback(16*Pp,1))
Ku=16,Pu=1.02-.639
Kp=0.6*Ku;Ki=1.2*Ku/Pu;Kd=0.075*Ku*Pu;C=Kp+Ki/s+Kd*s/(.01*s+1)
step(feedback(C*Pp,1))
bode(C*Pp)
grid