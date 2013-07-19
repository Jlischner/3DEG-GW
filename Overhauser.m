more off;
rs = 3.93;
kf = (9*pi/4)^(1/3)/rs; %# see Ashcroft (a.u.)
ks = [0.0: 0.1 :1.0]'*kf;
ef = kf^2/2;

eta = 0.15/27.21; %# broadening in chi0 in hartree
ws = [-5: 0.05 :2]'*ef; %# in hartrees
qs = [0.0001: 0.01 :10]'*kf;
us = [-1: 0.01 :1]';

Nfreq = length(ws);
Nq = length(qs);
Nk = length(ks);
Nu = length(us);
n = 3/(4*pi*rs^3);
wp =  sqrt(4*pi*n);
%# plasmon pole model from Lundqvist 2nd paper, Eq.7:
wq1 = sqrt( wp^2 + kf^2/3*qs.^2 + qs.^4/4);

x = qs/2/kf;
f = 0.5 + (1-x.^2)./(4*x).*log(abs( (1+x)./(1-x) ));
chi0 = -kf/pi^2 * f;
[beta,gamma,Gm_inf] = getGm_ov(rs);
Gm = (beta*x.^2 + gamma*x.^4)./(1+gamma*x.^4/Gm_inf);
vq = 4*pi./qs.^2;
Ixc = 1/2 * getIxc(n);
#Gm  = - Ixc./vq;

chiS = -chi0 ./ (1+Gm.*vq.*chi0);
wq_pa = sqrt(n./chiS) .* qs;

[alpha,Gp_inf] = getGp_ov(rs);
Gp = (1+alpha)*x.^2 ./(1 + (1+alpha)*x.^2/Gp_inf);
fxc = getfxc(n);
#Gp  = -fxc./vq;

epsO = 1 - vq.*chi0 ./(1 + Gp.*vq.*chi0);
wq = wp * sqrt( epsO ./ (epsO-1) );

Sig_PP   = zeros(Nfreq,Nk);
Sig_OV   = zeros(Nfreq,Nk); 
SigXs = zeros(Nk,1);

for ik = 1:Nk;
  printf("doing ik=%f out of Nk=%d \n",ik,Nk);
  k = ks(ik);
  
  for iw = 1:Nfreq;
    w = ws(iw);
    Sig_u = zeros(Nu,1);
    Sig_u2= zeros(Nu,1);

    for iu = 1:Nu;
      u = us(iu);
      xi = (k^2+qs.^2-2*k*qs*u)/2 - ef;
      Sig_u(iu) = trapz(qs,  (1-Gp*0).^2./(wq .* (w-xi-wq.*sign(xi) + I*eta)));    
      Sig_u2(iu)= trapz(qs, Gm.^2./(wq_pa .*(w-xi-wq_pa.*sign(xi) + I*eta)));
    endfor; % u-loop

      Sig_PP(iw,ik) += wp^2/(2*pi)*trapz(us, Sig_u );
      Sig_OV(iw,ik) += 6*n * trapz(us, Sig_u2);
  endfor; % w-loop

  %# bare exchange (Ashcroft Ch. 17, Eq. [17.19]):
  FX = 0.5 + (kf^2-k^2)/4/k/kf * log(abs( (kf+k)/(kf-k) ));
  if( abs(k-kf)  < 0.00001);
    FX = 0.5;
  elseif( abs(k) < 0.00001);
    FX = 1;
  endif;
  SigX = -2*kf/pi * FX;
  SigXs(ik) = SigX;
  Sig_PP(:,ik) += SigX;

endfor; %# ik-loop

%# get dmu:
[a,ikf] = min( abs( ks-kf ));
printf("using k=%f * kf as Fermi wave vector \n",ks(ikf)/kf);
dmu_PP = spline(ws,real(Sig_PP(:,ikf)),0);
dmu_OV = spline(ws,real(Sig_PP(:,ikf) + Sig_OV(:,ikf) ),0);

printf("rs = %f \n",rs);
printf("chemical potential mu = %f ry \n",dmu_PP*2);

DSig_PP = zeros(Nk,1);
DSig_OV = zeros(Nk,1);
for ik=1:Nk;
  ek = (ks(ik).^2/2 - ef);
  DSig_PP(ik) = interp1(ws,real(Sig_PP(:,ik)),ek);
  DSig_OV(ik) = interp1(ws,real(Sig_OV(:,ik)),ek);
endfor;

wss = [-5:0.01:1.5]'/27.21;
Nfreqs = length(wss);
A_PP = zeros(Nfreqs,Nk);
A_OV = zeros(Nfreqs,Nk);
for ii = 1:Nk;
  xi = ks(ii)^2/2 -ef;
  ReSig = interp1(ws,real(Sig_PP(:,ii)),wss);
  ImSig = interp1(ws,imag(Sig_PP(:,ii)),wss);
  A_PP(:,ii) = abs(ImSig)./( (wss-xi-ReSig+dmu_PP).^2 + ImSig.^2);

  ReSig = interp1(ws,real( Sig_PP(:,ii) + Sig_OV(:,ii) ),wss);
  ImSig = interp1(ws,imag(Sig_PP(:,ii) + Sig_OV(:,ii) ),wss);
  A_OV(:,ii) = abs(ImSig)./( (wss-xi-ReSig+dmu_OV).^2 + ImSig.^2);
endfor;
A_PP /= pi;
A_OV /= pi;

D2Sig_PP = zeros(Nk,1);
D2Sig_OV = zeros(Nk,1);
for ik=1:Nk;
  [a,b] = max(A_PP(:,ik) );
  D2Sig_PP(ik) = interp1(ws,real(Sig_PP(:,ik)),wss(b));
  
  [a,b] = max(A_OV(:,ik) );
  D2Sig_OV(ik) = interp1(ws,real(Sig_PP(:,ik)+Sig_OV(:,ik)),wss(b));
endfor;

printf("qp-energies: \n");
[a,b] = max(A_PP(:,1));
printf("PP: %f ry \n",wss(b)*2);
printf("H: %f \n",-ef*2);

more on;