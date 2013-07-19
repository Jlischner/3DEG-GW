eta = 0.001;
rs = 3.93;
kf = 1.92/rs; %# see Ashcroft (a.u.)
vf = kf; %# electrom mass = 1 in atomic units
ef = kf^2/2;
chiP = -kf/pi^2; %# static q=0 susceptibility for full system (=both spins)
n = 3/(4*pi*rs^3);
Ixc = getIxc(n); %# getIxc uses a.u.: computes Ixc for chi0/2
Ixc /= 2; %# get Ixc for full chi0 (=both spins)
fxc = getfxc(n);

qs = [0.01:0.01:40*kf]';
dq = qs(2)-qs(1);
Nq = length(qs);
ws = [0:0.001:1]'; %# in hartrees
Nfreq = length(ws);

chi0s = zeros(Nfreq,Nq);
chiSs = zeros(Nfreq,Nq);
chiCs = zeros(Nfreq,Nq);
chifs = zeros(Nfreq,Nq);
epsCs = zeros(Nfreq,Nq);
epsfs = zeros(Nfreq,Nq);
chiTs = zeros(Nfreq,Nq);

for iq = 1:Nq;

  q = qs(iq);

  [Ax,Ac,Bx,Bc]=Gtakada(rs);
  Gx = Ax*q.^2./(kf^2+Bx*q.^2);
  Gc = Ac*q.^2./(kf^2+Bc*q.^2);
  Gp = Gx + Gc;

  xi = q^2/2 - ef;
  x1 = q/kf + ( ws + I*eta)/ef/(q/kf);
  x2 = q/kf - ( ws + I*eta)/ef/(q/kf);
  
  f1 = (1 - x1.^2/4) .* log(  (x1+2)./(x1-2)  );
  f2 = (1 - x2.^2/4) .* log(  (x2+2)./(x2-2)  );
  
  u_qw = 0.5 + (f1 + f2)/(4*q/kf); %# u -> 1 (q ->0, w=0)
  chi0 = chiP * u_qw;
  chiS = chi0 ./(1-Ixc*chi0); %# interacting spin susceptibility

  vq = 4*pi/q^2;
  chiC = chi0 ./(1-vq*chi0); %# interacting charge susceptibility
  chif = chi0 ./(1-(vq+fxc)*chi0); %# adding fxc
  chit = chi0 ./(1-(1-Gp)*vq*chi0);

  epsC = 1+chiC*vq;
  epsf = 1+chif*vq;

  chi0s(:,iq) = chi0;
  chiSs(:,iq) = chiS;
  chiCs(:,iq) = chiC;
  chifs(:,iq) = chif;
  chiTs(:,iq) = chit;
  epsCs(:,iq) = epsC;
  epsfs(:,iq) = epsf;

end;

