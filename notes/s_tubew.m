function F=s_tubew(nw,n,N_w,N_p,F12,Fb,Fs);
%function F=s_tubew(nw,n,N_w,N_p,F12,Fb,Fs);
% computes view factor matrices for the furnace emulator

%KSTsakalis 11/02
%#inbounds
%#realonly
N_top=n-N_w-N_p; q1=min(nw-1,n); q3=nw+N_w+N_p-1;
F=[[F12(1:n,1:N_top);zeros(n,N_top)],zeros(2*n,2*n-N_top)];
if nw > N_top+1
F(N_top+1:q1,N_top+1:q1)=Fs(N_top+1:q1,N_top+1:q1);
end
F(N_top:q1,N_top)=F(N_top:q1,N_top)+Fb(1:q1-N_top+1);
if nw<n+2
F(1:q3,nw:q3)=F12(1:q3,nw:q3);
F(N_top:nw,nw)=F(N_top:nw,nw)+Fb([nw-N_top+1:-1:1]);
end    
