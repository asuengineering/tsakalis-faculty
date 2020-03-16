function [F1,F2,F3,F4,F5,w]=tube_vf(L,r,w,step_size);
% function [F1,F2,F3,F4,F5,w]=tube_vf(L,r,w,step_size);
% computation of view factors to tube slices
% F1 = inner to outer tube
% F2 = outer base to outer tube
% F3 = inner base to outer tube
% F4 = base to surrounding tube (F2|r=0)
% w = grid points

% KSTsakalis 11/02, 

if nargin <4; step_size=[]; end
if nargin <3; w=[]; end
if isempty(step_size);step_size=2e-3; end
step0=step_size; step1=step_size; 
if isempty(w); w=[0:1000]'/1000*L; end

D0=logspace(-3,log10(max(L,1)),40)';
k0=[0:step0:1]';  R0=[0:step1:1];
IN=0*k0+1;IN(1)=1/2;IN(length(IN))=1/2;
IN1=0*R0+1;IN1(1)=1/2;IN1(length(IN1))=1/2;

disp(' F1 computation: inner to outer tube')
D=w';
k=k0*acos(r);    K=k*(0*D+1);   DD=(0*k0+1)*D;
A=(1+r^2+DD.^2-2*r*cos(K)).^2;
%    A=max(1e-15,A);
F1=2*IN'*((r-cos(K)).*(r*cos(K)-1)./A)*k(2)/pi;
%vF1=sum(sum(FF1)*w(2))*w(2)/L;
%disp(['inner tube to outer tube  ',num2str(vF1)])

disp(' F5 computation: tube to self')
D=w';
k=k0*2*pi;    K=k*(0*D+1);   DD=(0*k0+1)*D;
A=(DD.^2+2*(1-cos(K))).^2;
   A=max(1e-15,A);
F5=2*IN'*((1-cos(K)).^2./A)*k(2)/2/pi;
F5(1)=1/2;

disp(' F2 computation: outer base to outer tube')
F2=0*D0;D=D0;
for i=1:length(D);
    C=-(2*D(i)/pi/(1-r^2));
    R=R0*(1-r)+r;
    K=k0*(acos(r)+acos(R));
    RR=(0*k0+1)*R;
    A=(1+RR.^2+D(i).^2-2*RR.*cos(K)).^2;
%    A=max(1e-15,A);
    F2(i)=C*(IN'*(RR.*(RR.*cos(K)-1)./A)).*K(2,:)*IN1'*(R(2)-R(1));
end
F2 = INTERP1(D,2*F2,w,'linear','extrap');
%vF2=sum(F2)*w(2)*(1-r^2)/2/L;
%disp(['outer tube to outer base   ',num2str(vF2)])
%disp(['outer tube to self   ',num2str(1-vF1*r-2*vF2)])

disp(' F3 computation: outer base to inner tube')
F3=0*D0;D=D0;
for i=1:length(D);
    C=(2*D(i)*r/pi/(1-r^2));
    R=R0*(1-r)+r;
    K=k0*(acos(r./R));
    RR=(0*k0+1)*R;
    A=(r^2+RR.^2+D(i).^2-2*r*RR.*cos(K)).^2;
%    A=max(1e-15,A);
    F3(i)=C*(IN'*(RR.*(RR.*cos(K)-r)./A)).*K(2,:)*IN1'*(R(2)-R(1));
end
F3 = INTERP1(D,2*F3,w,'linear','extrap');
%vF3=sum(F3)*w(2)*(1-r^2)/2/r/L;
%disp(['inner tube to outer base   ',num2str(vF3)])

disp(' F4 computation: entire base to surrounding tube') 
F4=0*D;D=D0;r=0;
for i=1:length(D);
    R=R0*(1-r)+r;
    K=k0*(acos(r)+acos(R));
    RR=(0*k0+1)*R;
    A=(1+RR.^2+D(i).^2-2*RR.*cos(K)).^2;
%    A=max(1e-15,A);
    F4(i)=-(2*D(i)/pi/(1-r^2))*(IN'*(RR.*(RR.*cos(K)-1)./A)).*K(2,:)*IN1'*(R(2)-R(1));
end
F4 = INTERP1(D,2*F4,w,'linear','extrap');
%vF4=sum(F4)*w(2)*(1-0^2)/L/2;
%disp(['outer tube to entire base   ',num2str(vF4)])


