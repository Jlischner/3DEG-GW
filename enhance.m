rss = [0.1: 0.1 : 6.0]';
Nrs = length(rss);
enh = zeros(Nrs,1);
enh2= zeros(Nrs,1);
enh3= zeros(Nrs,1);
kappa = zeros(Nrs,1);
g = zeros(Nrs,1);

for irs = 1:Nrs;

  rs = rss(irs);
  kf = 1.92/rs; %# see Ashcroft (a.u.)
  chiP = -kf/pi^2; %# static q=0 susceptibility for full system (=both spins)
  n = 3/(4*pi*rs^3);
  Ixc = getIxc(n); %# getIxc uses a.u.: computes Ixc for chi0/2
  Ixc /= 2; %# get Ixc for full chi0 (=both spins)
  printf("rs=%f I_lda=%f \n",rs,Ixc);
  enh(irs) = 1./(1-Ixc*chiP); %# interacting spin susceptibility

  [Ax,Ac,Bx,Bc,kappaF,g0]=Gtakada(rs);
  #enh3(irs) = enh_tak;
  kappa(irs) = kappaF;
  g(irs) = g0;
  Ixc = 4*pi*(Ax-Ac)/kf^2;
  printf("rs=%f I_tak=%f \n",rs,Ixc);
  enh2(irs) = 1./(1-Ixc*chiP);
end;

