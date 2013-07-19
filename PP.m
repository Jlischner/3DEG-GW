rs = 3.93;
kf = (9*pi/4)^(1/3)/rs; %# see Ashcroft (a.u.)
ks = [0.0: 0.1 :1.0]'*kf;
ef = kf^2/2;

eta = 0.15/27.21; %# broadening in chi0 in hartree
ws = [-20: 0.05 :20]'*ef; %# in hartrees
qs = [0.0001: 0.01 :10]'*kf;
us = [-1: 0.01 :1]';

Nfreq = length(ws);
Nq = length(qs);
Nk = length(ks);
Nu = length(us);
n = 3/(4*pi*rs^3);
wp =  sqrt(4*pi*n);
%# plasmon pole model from Lundqvist 2nd paper, Eq.7:
wq = sqrt( wp^2 + kf^2/3*qs.^2 + qs.^4/4);

Sig_PP   = zeros(Nfreq,Nk);
ReSig_PP = zeros(Nfreq,Nk);
SigXs = zeros(Nk,1);

for ik = 1:Nk;
  k = ks(ik);
  
  for iw = 1:Nfreq;
    w = ws(iw);
    Sig_u = zeros(Nu,1);
    
    for iu = 1:Nu;
      u = us(iu);
      xi = (k^2+qs.^2-2*k*qs*u)/2 - ef;
      Sig_u(iu) = trapz(qs, 1./(wq .* (w-xi-wq.*sign(xi) + I*eta)));    
    endfor; % u-loop

      Sig_PP(iw,ik) += wp^2/(2*pi)*trapz(us, Sig_u );
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

  etaKK = 0.001/27.21;
  ReSig_PP(:,ik) = getReSig(ws,-imag(Sig_PP(:,ik)),etaKK);
  ReSig_PP(:,ik) += SigX;

endfor; %# ik-loop

%# get dmu:
[a,ikf] = min( abs( ks-kf ));
printf("using k=%f * kf as Fermi wave vector \n",ks(ikf)/kf);
dmu_PP = spline(ws,real(Sig_PP(:,ikf)),0);
dmu_PP2 = spline(ws,ReSig_PP(:,ikf),0);

printf("rs = %f \n",rs);
printf("chemical potential mu = %f ry \n",dmu_PP*2);
printf("chemical potential mu2= %f ry \n",dmu_PP2*2);

DSig_PP = zeros(Nk,1);
for ik = 1:Nk;
  ek = (ks(ik).^2/2 - ef);
  DSig_PP(ik) = interp1(ws,real(Sig_PP(:,ik)),ek);
endfor;

wss = [-8:0.01:0.5]'/27.21;
Nfreqs = length(wss);
A_PP = zeros(Nfreqs,Nk);
A_PP2= zeros(Nfreqs,Nk);
for ii = 1:Nk;
  xi = ks(ii)^2/2 -ef;
  ReSig = interp1(ws,real(Sig_PP(:,ii)),wss);
  ImSig = interp1(ws,imag(Sig_PP(:,ii)),wss);
  A_PP(:,ii) = abs(ImSig)./( (wss-xi-ReSig+dmu_PP).^2 + ImSig.^2);

  ReSig = interp1(ws,ReSig_PP(:,ii),wss);
  A_PP2(:,ii) = abs(ImSig)./( (wss-xi-ReSig+dmu_PP2).^2 + ImSig.^2);
endfor;
A_PP /= pi;
A_PP2 /= pi;

D2Sig_PP = zeros(Nk,1);
for ik = 1:Nk;
  [a,b] = max(A_PP(:,ik) );
  D2Sig_PP(ik) = interp1(ws,real(Sig_PP(:,ik)),wss(b));
endfor;

printf("qp-energies: \n");
[a,b] = max(A_PP(:,1));
printf("PP: %f ry \n",wss(b)*2);
[a,b] = max(A_PP2(:,1));
printf("PP2: %f ry \n",wss(b)*2);
printf("H: %f \n",-ef*2);