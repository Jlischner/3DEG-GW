rs = 4.00;
kf = 1.92/rs; %# see Ashcroft (a.u.)
ks = [0.0: 1.0 :1.0]'*kf;

eta = 0.001; %# broadening in chi0 in hartree
ws = [-2: 0.01 :2]'; %# in hartrees
qs = [0.01: 0.01 :5*kf]';
us = [-1: 0.01 :1]';

Nfreq = length(ws);
dq = qs(2)-qs(1);
Nq = length(qs);
Nk = length(ks);
du = us(2)-us(1);
Nu = length(us);
vf = kf; %# electron mass = 1 in atomic units
ef = kf^2/2;
chiP = -kf/pi^2; %# static q=0 susceptibility for full system (=both spins)
n = 3/(4*pi*rs^3);
Ixc = getIxc(n); %# getIxc uses a.u.: computes Ixc for chi0/2
Ixc /= 2; %# get Ixc for full chi0 (=both spins)
fxc = getfxc(n);
wp =  sqrt(4*pi*n);

ImSig_SF = zeros(Nfreq,Nk);
ImSig_RP = zeros(Nfreq,Nk);
ImSig_OK = zeros(Nfreq,Nk);
ImSig_LH = zeros(Nfreq,Nk);

ReSig_SF = zeros(Nfreq,Nk);
ReSig_RP = zeros(Nfreq,Nk);
ReSig_OK = zeros(Nfreq,Nk);
ReSig_LH = zeros(Nfreq,Nk);
ReSig_PP = zeros(Nfreq,Nk);
Sig_PP   = zeros(Nfreq,Nk);
SigXs = zeros(Nk,1);

for ik = 1:Nk;
  k = ks(ik);
  
  for iq = 1:Nq;
    for iu = 1:Nu;
      
      q = qs(iq);
      u = us(iu);
      xi = (k^2+q^2-2*k*q*u)/2 - ef;
      wshift = sign(xi)*(ws-xi);

      x1 = q/kf + ( wshift + I*eta)/ef/(q/kf);
      x2 = q/kf - ( wshift + I*eta)/ef/(q/kf);
      
      f1 = (1 - x1.^2/4) .* log(  (x1+2)./(x1-2)  );
      f2 = (1 - x2.^2/4) .* log(  (x2+2)./(x2-2)  );
      
      u_qw = 0.5 + (f1 + f2)/(4*q/kf); %# u -> 1 (q ->0, w=0)
      chi0 = chiP * u_qw; %# im(chi0) is anti-symmetric
      chiS = chi0 ./(1-Ixc*chi0); %# interacting spin susceptibility
      
      vq = 4*pi/q^2;
      vf = vq + fxc;
      chiC = chi0 ./(1-vq*chi0); %# interacting RPA charge susceptibility
      chif = chi0 ./(1-vf*chi0); %# interacting LDA charge susceptibility
      
      ImSig_SF(:,ik) += Ixc^2*q^2 * imag(chiS) .* ( wshift > 0);
      ImSig_RP(:,ik) += vq^2 *q^2 * imag(chiC) .* ( wshift > 0);
      ImSig_OK(:,ik) += vf^2 *q^2 * imag(chif) .* ( wshift > 0);
      ImSig_LH(:,ik) += vq^2 *q^2 * imag(chif) .* ( wshift > 0);
      
      wq = ef*sqrt( wp^2/ef^2 + 16/3*q^2/kf^2 + q^4/kf^4);
      Sig_PP(:,ik) += q^2 * vq * wp^2/wq * 1./(ws-xi-wq*sign(xi) + I*eta);

    endfor;
  endfor;

  ImSig_SF(:,ik) *= -3/(2*pi)^2*dq*du;
  ImSig_RP(:,ik) *= -1/(2*pi)^2*dq*du;
  ImSig_OK(:,ik) *= -1/(2*pi)^2*dq*du;
  ImSig_LH(:,ik) *= -1/(2*pi)^2*dq*du;
  Sig_PP(:,ik)   *=  1/(2*pi)^2*dq*du;

  %# bare exchange (Ashcroft Ch. 17, Eq. [17.19]):
  FX = 0.5 + (kf^2-k^2)/4/k/kf * log(abs( (kf+k)/(kf-k) ));
  if( abs(k-kf)  < 0.00001);
    FX = 0.5;
  elseif( abs(k) < 0.00001);
    FX = 1;
  endif;
  SigX = -2*kf/pi * FX;
  SigXs(ik) = SigX;

  etaKK = 0.001/27.21;
  ReSig_SF(:,ik) = getReSig(ws,ImSig_SF(:,ik),etaKK);
  ReSig_RP(:,ik) = getReSig(ws,ImSig_RP(:,ik),etaKK) + SigX;
  ReSig_OK(:,ik) = getReSig(ws,ImSig_OK(:,ik),etaKK) + SigX;
  ReSig_LH(:,ik) = getReSig(ws,ImSig_LH(:,ik),etaKK) + SigX;
  ReSig_PP(:,ik) = getReSig(ws,-imag(Sig_PP(:,ik)),etaKK) + SigX;
  Sig_PP(:,ik) += SigX;

endfor; %# ik-loop

ws *= 27.21;
wp *= 27.21;
ImSig_SF *= 27.21;
ImSig_RP *= 27.21;
ImSig_OK *= 27.21;
ImSig_LH *= 27.21;
ReSig_SF *= 27.21;
ReSig_RP *= 27.21;
ReSig_OK *= 27.21;
ReSig_LH *= 27.21;
Sig_PP   *= 27.21;
ReSig_PP *= 27.21;

%# get dmu:
[a,ikf] = min( abs( ks-kf ));
printf("using k=%f * kf as Fermi wave vector \n",ks(ikf)/kf);
dmu_RP = spline(ws,ReSig_RP(:,ikf),0);
dmu_LH = spline(ws,ReSig_LH(:,ikf),0);
dmu_PP = spline(ws,real(Sig_PP(:,ikf)),0);

%# calculate on-shell self energy:
DSig_OK = zeros(Nk,1);
DSig_RP = zeros(Nk,1);
DSig_SF = zeros(Nk,1);
for ii = 1:Nk;
  xi = ( ks(ii)^2/2 -ef )*27.21;
  DSig_OK(ii) = spline(ws,ReSig_OK(:,ii),xi);
  DSig_RP(ii) = spline(ws,ReSig_RP(:,ii),xi);
  DSig_SF(ii) = spline(ws,ReSig_SF(:,ii),xi);
endfor;

wss = [-20:0.01:20]';
Nfreqs = length(wss);
A_RP = zeros(Nfreqs,Nk);
A_LH = zeros(Nfreqs,Nk);
A_PP = zeros(Nfreqs,Nk);
for ii = 1:Nk;
  xi = (ks(ii)^2/2 -ef)*27.21;
  ReSig = interp1(ws,ReSig_RP(:,ii),wss);
  ImSig = interp1(ws,ImSig_RP(:,ii),wss);
  A_RP(:,ii) = abs(ImSig)./( (wss-xi-ReSig+dmu_RP).^2 + ImSig.^2);

  ReSig = interp1(ws,ReSig_LH(:,ii),wss);
  ImSig = interp1(ws,ImSig_LH(:,ii),wss);
  A_LH(:,ii) = abs(ImSig)./( (wss-xi-ReSig+dmu_LH).^2 + ImSig.^2);

  ReSig = interp1(ws,real(Sig_PP(:,ii)),wss);
  ImSig = interp1(ws,imag(Sig_PP(:,ii)),wss);
  A_PP(:,ii) = abs(ImSig)./( (wss-xi-ReSig+dmu_PP).^2 + ImSig.^2);
endfor;
A_RP /= pi;
A_LH /= pi;
A_PP /= pi;

plot(ws/wp,ReSig_RP(:,1)/wp,'sr-',ws/wp,ImSig_RP(:,1)/wp,'sb-');axis([-3 3 -5 3]);