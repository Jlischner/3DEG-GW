eta = 0.001;
rs = 4.0;
kf = 1.92/rs; %# see Ashcroft (a.u.)
vf = kf; %# electrom mass = 1 in atomic units
ef = kf^2/2;
chiP = -kf/pi^2; %# static q=0 susceptibility for full system (=both spins)
n = 3/(4*pi*rs^3);
Ixc = getIxc(n); %# getIxc uses a.u.: computes Ixc for chi0/2
Ixc /= 2; %# get Ixc for full chi0 (=both spins)

qs = [0.01:0.01:10*kf]';
dq = qs(2)-qs(1);
Nq = length(qs);
ws = [-1:0.02:1]'; %# in hartrees
Nfreq = length(ws);

ImSig_SF = zeros(Nfreq,1);
ImSig_GW = zeros(Nfreq,1);

for iq = 1:Nq;

  q = qs(iq);
  xi = q^2/2 - ef;
  wshift = sign(xi)*(ws-xi);

  x1 = q/kf + ( wshift + I*eta)/ef/(q/kf);
  x2 = q/kf - ( wshift + I*eta)/ef/(q/kf);
  
  f1 = (1 - x1.^2/4) .* log(  (x1+2)./(x1-2)  );
  f2 = (1 - x2.^2/4) .* log(  (x2+2)./(x2-2)  );
  
  u_qw = 0.5 + (f1 + f2)/(4*q/kf); %# u -> 1 (q ->0, w=0)
  chi0 = chiP * u_qw; %# im(chi0) is anti-symmetric
  chiS = chi0 ./(1-Ixc*chi0); %# interacting spin susceptibility

  vq = 4*pi/q^2;
  chiC = chi0 ./(1-vq*chi0); %# interacting charge susceptibility
  
  ImSig_SF += Ixc^2*q^2    * imag(chiS) .* ( wshift > 0);
  ImSig_GW += (4*pi)^2/q^2 * imag(chiC) .* ( wshift > 0);

end;

ImSig_SF = -3/(2*pi^2)*dq * ImSig_SF;
ImSig_GW = -1/(2*pi^2)*dq * ImSig_GW;
ws *= 27.21;
ImSig_SF *= 27.21;
ImSig_GW *= 27.21;

eta = 0.01; %# in eV
ReSig_SF = getReSig(ws,ImSig_SF,eta);
plot(ws,ImSig_SF,'g-',ws,ImSig_GW,'r-',ws,ReSig_SF,'m-')