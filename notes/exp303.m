disp('Experiment 1')
%Exp.1: Signal generation, Fourier Transform
fs=44100; %sample sin at 44.1 kHz
no_pts=2*8192;
t=([1:no_pts]'-1)/fs;
y1=sin(2*pi*1000*t);
plot(t,y1);axis([0,0.004,-1.2,1.2]);title('The signal y(t)');pause
disp('original');sound(y1,fs);pause(no_pts/fs*2)
fr=([1:no_pts]'-1)/no_pts*fs;fr=fr(1:no_pts/2); %in Hz
f1=abs(fft(y1));f1=f1(1:no_pts/2)/fs;
frp=fr*2*pi;tmax=max(t);
F1=1/j*sin((frp-1000*2*pi)*tmax/2).*exp(j*(frp-1000*2*pi)*tmax/2)./(frp-1000*2*pi);
F2=1/j*sin((frp+1000*2*pi)*tmax/2).*exp(j*(frp+1000*2*pi)*tmax/2)./(frp+1000*2*pi);
F=abs(F1-F2);
loglog(fr,F,fr,f1);title('Frequency Domain magnitudes: F[y], FFT[y]');pause

disp('Experiment 2')
%Exp.2: Low-frequency sampling and reconstruction, Fourier Transform
a=4; %sample sin at 44.1/a kHz; use a = power of 2 
t_a=([1:no_pts/a]'-1)/fs*a;
y_a=sin(2*pi*1000*t_a);
plot(t,y1,t_a,y_a);axis([0,0.004,-1.2,1.2]);
title('y(t) sampled at high and low rates');pause
disp('original');sound(y1,fs);pause(no_pts/fs*2)
disp('low-freq sampling');sound(y_a,fs/a);pause(no_pts/fs*2)
% this sounds OK because of the internal filtering
fr_a=([1:no_pts/a]'-1)/no_pts*a*fs/a;fr_a=fr_a(1:no_pts/a/2);
f_a=abs(fft(y_a));f_a=f_a(1:no_pts/a/2)/fs*a;
loglog(fr,F,fr,f1,fr_a,f_a);title('Magnitudes of F[y], FFT[y], FFT[y_a]');pause

%bypass internal filters and emulate 44.1/a kHz D/A
y_n = interp1(t_a,y_a,t,'nearest','extrap');
        % REM: nearest neighbor interpolation ~ shifted ZOH
f_n=abs(fft(y_n));f_n=f_n(1:no_pts/2)/fs;
plot(t,y1,t,y_n);axis([0,0.004,-1.2,1.2])
title('Original and low-rate ZOH signals');pause
disp('original');sound(y1,fs);pause(no_pts/fs*2)
disp('nearest (~ZOH)');sound(y_n,fs);pause(no_pts/fs*2)
loglog(fr,F,fr,f1,fr,f_n);title('Magnitudes of F[y], FFT[y], FFT[y_{a,ZOH}]');pause

disp('Experiment 3')
%Exp.3: Digital filtering and reconstruction from the low-rate sampled signal
%linear interpolation
y_l = interp1(t_a,y_a,t,'linear','extrap');
f_l=abs(fft(y_l));f_l=f_l(1:no_pts/2)/fs;
plot(t,y1,t,y_l);axis([0,0.004,-1.2,1.2])
title('Original and lin.-interpolated low-rate signals');pause
disp('original');sound(y1,fs);pause(no_pts/fs*2)
disp('lin.interp.');sound(y_n,fs);pause(no_pts/fs*2)
loglog(fr,F,fr,f1,fr,f_l);title('Magnitudes of F[y], FFT[y], FFT[y_l]');pause

%upsampling
y_u=y1*0;k=1:a:no_pts;y_u(k)=y_a;
f_u=abs(fft(y_u));f_u=f_u(1:no_pts/2)/fs*a;
plot(t,y1,t,y_u);axis([0,0.004,-1.2,1.2])
title('Original and upsampled signals');pause
disp('original');sound(y1,fs);pause(no_pts/fs*2)
disp('upsampled');sound(y_u,fs);pause(no_pts/fs*2)
loglog(fr,F,fr,f1,fr,f_u);title('Magnitudes of F[y], FFT[y], FFT[y_u]');pause

%3kHz lowpass filter (1000*2-1 point approximation)
hs=sin(3000*2*pi*t)/pi./t; % this gives a division by zero warning
i=1000:-1:2;HS=[hs(i);3000*2;hs(2:1000)]*a/fs;  % set the correct value at 0 
y_f=conv(HS,y_u);
%eliminate initial/final points (initialization) of 1000 samples = 2.27e-2 s (significant)
y_f=y_f(1000:999+no_pts);
f_f=abs(fft(y_f));f_f=f_f(1:no_pts/2)/fs;
plot(t,y1,t,y_f);axis([0,0.004,-1.2,1.2])
title('Original and filtered upsampled signals');pause
disp('original');sound(y1,fs);pause(no_pts/fs*2)
disp('filtered-upsampled');sound(y_f,fs);pause(no_pts/fs*2)
loglog(fr,F,fr,f1,fr,f_f);title('Magnitudes of F[y], FFT[y], FFT[y_f]');pause

%use a butterworth filter
A=[ 1.0000e+000 -7.5419e+000  2.5836e+001 -5.2896e+001  7.1625e+001 -6.6986e+001 ...  
    4.3797e+001 -1.9759e+001  5.8842e+000 -1.0441e+000  8.3813e-002];
B=[  2.2721e-008  2.2721e-007  1.0225e-006  2.7265e-006  4.7714e-006  5.7257e-006 ...
    4.7714e-006  2.7265e-006  1.0225e-006  2.2721e-007  2.2721e-008];
td=-3.7885e-4;
% comment the following two calls if the functions are not available
[B,A] = butter(10,2700/22050);
[m,p]=dbode(B,A,1/fs,2*pi*1000); td=p/180/2/1000; % find delay at 1kHz
y_b=dlsim(B*a,A,y_u);
f_b=abs(fft(y_b));f_b=f_b(1:no_pts/2)/fs;
plot(t,y1,t+td,y_b);axis([0,0.004,-1.2,1.2]);
title('Original and butterworth-filtered upsampled signals');pause
disp('original');sound(y1,fs);pause(no_pts/fs*2)
disp('b-filtered upsampled');sound(y_b,fs);pause(no_pts/fs*2)
loglog(fr,F,fr,f1,fr,f_b);title('Magnitudes of F[y], FFT[y], FFT[y_b]');pause



%Exp.4: Equalization
[y,f]= wavread('tada');
y=y(:,1); % keep only one channel
ly=length(y); ty=([1:ly]'-1)/f;
w=([1:ly]'-1)/ly*f;w=w(1:ly/2); %in Hz
fy=abs(fft(y));fy=fy(1:ly/2)/f; % FFT
    soundsc(y,f)
loglog(w,fy);title('Magnitude of F[y]');pause

i=2:500;i_=500:-1:2; L=500:499+ly; % define 3 1000-point lowpass filters 
hs1=sin(700*2*pi*ty)/pi./ty;HS1=[hs1(i_);700*2;hs1(i)]/f;
hs2=sin(2000*2*pi*ty)/pi./ty;HS2=[hs2(i_);2000*2;hs2(i)]/f;
hs3=sin(5000*2*pi*ty)/pi./ty;HS3=[hs3(i_);5000*2;hs3(i)]/f;
z1=conv(HS1,y);z1=z1(L);
z2=conv(HS2,y);z2=z2(L);
z3=conv(HS3,y);z3=z3(L);
done=0;
while done ==0
    coef=input('enter equalization coefficients [low, mid-low, mid-hi, hi] ');
    if isempty(coef),coef=[1 1 1 1];end
    yx=[z1 z2-z1 z3-z2 y-z3]*coef';
    soundsc(yx,f)
    fx=abs(fft(yx));fx=fx(1:ly/2)/f;
    loglog(w,fx);title('Equalized signal');pause
    done=input('Done 1=y, [0]  ');
    if isempty(done);done=0;end
end

