% Initialization of the wafer furnace emulator
% Use lower sample times in case of wafer/gas temperature instability
% primary heat loss adjustment through lam_bc and f_amb parameters
% other adjustables: pow_pml, DX_T, DX_P, top, ped mass/conductivity adjustments
%     n_sl: 2 wafers/slice ~5x real-time, 5-10 deg errors (but similar behavior)
%     n_sl: 1, ~5/6 real-time
%     without direct element-wafer radiation, ~2x faster.

% KSTsakalis 11/02, v1.2

n_sl=2;                      % number of wafers per slice
N_w=100/n_sl;                % number of wafers  per slice 
D_w=0.2;                     % Wafer diameter
dw_t = 0.001*n_sl;           % wafer thickness, add some for boat mass
dx_w=0.01*n_sl;              % intra-wafer distance
DX_T = 0.15;                 % top space (shower head)
DX_P = 0.20;                 % pedestal height
DX_R = 0.04;                 % wafer-tube distance
DX_E = 0.05;                 % element-tube distance
T_01=500;T_a01=25;T_amb=25;   % initial and ambient temperatures 
d_ = 0.007;                  % tube thickness
d_e = 0.006;                 % element thickness
d_g = 0.006;                  % element gaps
d_ins = 0.2;                 % outer insulation thickness
elem_tmax = 1250;
tube_tmax = 1200;
pow_pml=200;                 % power per unit length

Tsa=.5; decim=10;             % sampling time, data collection decimation

Li = dx_w;                   % length resolution

N_top = round(DX_T/Li);  N_ped = round(DX_P/Li); 
DX_T = N_top*Li;  DX_P = N_ped*Li; n = N_w+N_top+N_ped;  L = Li;

T_0=T_01*ones(n,1); T_a0=T_a01*ones(n,1); T_w0 = T_amb*ones(N_w,1);
T_w0_top = T_01*ones(N_top,1);  T_w0_ped = T_amb*ones(N_ped,1);

rho_e = 8500;  e_e = 0.8;  rho_p = 2400;  e_p = 0.9; e_ew=0.1;
rho_w = 2330;  e_w = 0.6;  pi=3.14159; h_gas = 2;

D_e = D_w+2*DX_R+2*d_+2*DX_E;
D_po = D_w+2*DX_R+2*d_;  
D_pi = D_w+2*DX_R;
A_e = pi*D_e*L;  V_e = pi*D_e*L*d_e/(d_e+d_g)*d_e; m_e = rho_e*V_e;
pow_pm = pow_pml*pi*D_e*n*L/(d_e+d_g);
A_amb = pi*(D_e+2*d_ins)*L;  h_amb = 1;
V_s = pi/4*(D_e^2-D_po^2)*L;  
V_p = pi/4*(D_po^2-D_pi^2)*L;  m_p = V_p*rho_p;  
A_pi = pi*D_pi*L;  
A_po = pi*D_po*L;
V_a = pi/4*D_pi^2*L;
A_w = pi/4*D_w^2*n_sl; V_w = pi/4*D_w^2*dw_t;
m_w = rho_w*V_w;
A_wc=pi*D_w*Li;

mcp_adj_top = 0.25;  k_adj_top = 0.5;
V_top = pi*D_w^2/4*Li;  m_top = rho_p*V_top;
f_amb_top = .3; 
A_amb_top = pi*D_w^2/4; 

mcp_adj_ped = 0.36;  k_adj_ped = 0.6;
V_ped = pi*D_w^2/4*Li;  m_ped = V_ped*rho_p; 
f_amb_ped = 1;
A_amb_ped = pi*D_w^2/4; 

lam_p_bct=.4; lam_p_bcb=1; 

SW = 1+(4*dx_w^2+D_w^2)/(D_w^2); FW = (SW-sqrt(SW^2-4))/2;
S_t = 1+(4*(Li)^2+D_pi^2)/(D_w^2); F_t = (S_t-sqrt(S_t^2-4*(D_pi/D_w)^2))/2;
S_b = 1+(4*(Li)^2+D_pi^2)/(D_w^2); F_b = (S_b-sqrt(S_b^2-4*(D_pi/D_w)^2))/2;
gww_bot = ones(N_w,1)*SW; gww_top = gww_bot;
gww_bot(N_w) = F_b; gww_top(1) = F_t;

%Element radiation 
len=n*Li;  ro=D_e/2;  ri=D_po/2;  ei=e_p; eo=e_e; w=[0:n-1]'/n; w0=w;
[F1,F2,F3,F4,F5]=tube_vf(len/ro,ri/ro,w*len,1e-2);
A_ep=pi*D_po*Li;
F_A_pe=A_ep./((1/ei-1)+(1/eo-1)*D_po/D_e+1./(toeplitz(F1*Li)));

Aob = pi*(D_e^2-D_po^2)/4;  Aib = pi*D_po^2/4;  A = pi*D_e*Li;
F_A_e = Aob./((1/ei-1)+(1/eo-1)*Aob/A+1./(F2*Li));
F_A_e=lam_p_bct*F_A_e+lam_p_bcb*F_A_e([n:-1:1]);

%Profile radiation
A=pi*D_po*Li;
F_A_p = Aob./((1/ei-1)+(1/ei-1)*Aob/A+1./(F3*Li));
F_A_p=lam_p_bct*F_A_p+lam_p_bcb*F_A_p([n:-1:1]);

%Profile-wafer radiation
len=2*n*Li;  ro=D_pi/2;  ri=D_w/2;  ei=e_w; eo=e_p; w=[0:2*n-1]'/n/2;
[F1,F2,F3,F4,F5]=tube_vf(len/ro,ri/ro,w*len,1e-2);
A=pi*D_pi*Li;  Aib=pi*D_w^2/4;  Aob = pi*(D_pi^2-D_w^2)/4; 
Fb_A_pw = Aib./((1/ei-1)+(1/eo-1)*Aib/A+1./(F4*Li));  % approx inner base

A2=pi*D_w*Li;
F_A_p2 = Aob./((1/ei-1)+(1/ei-1)*Aob/A2+1./(F2*Li));
%F_A_p=F_A_p+F_A_p2(1:n)+F_A_p2([n:-1:1]);

A_pw=pi*D_w*Li;
F_A_pw = A_pw./((1/ei-1)+(1/eo-1)*D_pi/D_w+1./(toeplitz(F1*Li)));
Fs_p = A_pi./((1/eo-1)+(1/eo-1)+1./(toeplitz(F5*Li)));

%Element-wafer direct radiation
len=2*n*Li;  ro=D_e/2;  ri=D_w/2;  ei=e_w; eo=e_e; w=[0:2*n-1]'/n/2;
[F1,F2,F3,F4,F5]=tube_vf(len/ro,ri/ro,w*len,1e-2);
A=pi*D_e*Li;  Aib=pi*D_w^2/4;  Aob = pi*(D_e^2-D_w^2)/4; 
Fb_A_ew = e_ew*Aib./((1/ei-1)+(1/eo-1)*Aib/A+1./(F4*Li));  % approx inner base

A2=pi*D_w*Li;
F_A_e2 = Aob./((1/ei-1)+(1/ei-1)*Aob/A2+1./(F2*Li));
%F_A_p=F_A_p+F_A_p2(1:n)+F_A_p2([n:-1:1]);

A_ew=pi*D_w*Li;
F_A_ew = e_ew*A_ew./((1/ei-1)+(1/eo-1)*D_pi/D_w+1./(toeplitz(F1*Li)));
Fs_e = e_ew*A_e./((1/eo-1)+(1/eo-1)+1./(toeplitz(F5*Li)));

% TC measurement computations
n24=fix(N_w/3);
h_zones = [N_top,n24,N_w-2*n24,n24,N_ped];
spike_loc = [round(h_zones(1)/2),h_zones(1)+round(h_zones(2)/2),...
  sum(h_zones(1:2))+round(h_zones(3)/2),sum(h_zones(1:3))+round(h_zones(4)/2),...
  sum(h_zones(1:4))+round(h_zones(5)/2)];
profile_loc = [N_top,h_zones(1)+round(h_zones(2)/2),...
  sum(h_zones(1:2))+round(h_zones(3)/2),sum(h_zones(1:3))+round(h_zones(4)/2),...
  n-N_ped];

sp_e=0.3*DX_E; pr_t=0.4*DX_R; dw=[0:n-1]'*Li;
sp_t=DX_E-sp_e; pr_w=DX_R-pr_t;  ov=0*dw+1;

SPE=toeplitz(sp_e./(sqrt(sp_e^2+dw.^2).^3));
SPT=toeplitz(sp_t./(sqrt(sp_t^2+dw.^2).^3));
PRT=toeplitz(pr_t./(sqrt(pr_t^2+dw.^2).^3));
PRW=toeplitz(pr_w./(sqrt(pr_w^2+dw.^2).^3));
SPE=SPE(:,spike_loc);
SPT=SPT(:,spike_loc);  spsum=sum(SPE)+sum(SPT);
PRT=PRT(:,profile_loc);
PRW=PRW(:,profile_loc);  prsum=sum(PRT)+sum(PRW);

SPE=SPE.*(ov*(1./spsum));SPE=SPE';
SPT=SPT.*(ov*(1./spsum));SPT=SPT';
PRT=PRT.*(ov*(1./prsum));PRT=PRT';
PRW=PRW.*(ov*(1./prsum));PRW=PRW';

pr_gas = 0.2;
PRG=0*SPE; for ii=1:length(profile_loc); PRG(ii,profile_loc(ii))=1;end
H_Z = zeros(n,5);
H_Z(1:h_zones(1),1)=ones(h_zones(1),1);
H_Z(1+h_zones(1):sum(h_zones([1:2])),2)=ones(h_zones(2),1);
H_Z(1+sum(h_zones([1:2])):sum(h_zones([1:3])),3)=ones(h_zones(3),1);
H_Z(1+sum(h_zones([1:3])):sum(h_zones([1:4])),4)=ones(h_zones(4),1);
H_Z(1+sum(h_zones([1:4])):sum(h_zones),5)=ones(h_zones(5),1);

