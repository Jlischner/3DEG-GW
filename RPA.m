more off;

rs = 3.93;
kf = (9*pi/4)^(1/3)/rs; %# see Ashcroft (a.u.)
ks = [0.0: 0.1 :1.0]'*kf;
ef = kf^2/2;

eta = 0.15/2/27.21; %# broadening in chi0 in hartree
ws = [-15: 0.05 :30]'*ef; %# in hartrees
qs = [0.0001: 0.01 :10]'*kf;
us = [-1: 0.01 :1]';

Nfreq = length(ws);
Nq = length(qs);
Nk = length(ks);
Nu = length(us);
n = 3/(4*pi*rs^3);
wp =  sqrt(4*pi*n);
chiP = -kf/pi^2;
%# plasmon pole model from Lundqvist 2nd paper, Eq.7:
wq = sqrt( wp^2 + kf^2/3*qs.^2 + qs.^4/4);
vq = 4*pi./qs.^2;

%# LDA local field factors:
#printf("using LDA local field factors \n");
#Ixc = getIxc(n); %# getIxc uses a.u.: computes Ixc for chi0/2
#Ixc /= 2; %# get Ixc for full chi0 (= both spins)
#fxc = getfxc(n);
#vf = vq + fxc;
%#*************************

%# Takada local field factors:
#printf("using Takada local field factors \n");
#[Ax,Ac,Bx,Bc]=Gtakada(rs);
#Gx = Ax*qs.^2./(kf^2+Bx*qs.^2);
#Gc = Ac*qs.^2./(kf^2+Bc*qs.^2);
#Gp = Gx + Gc;
#Gm = Gx - Gc;
#vf = (1-Gp).* vq;
#Ixc= Gm .* vq;
%#**************************

%# Overhauser local field factors:
printf("Using Overhauser local field factors \n");
x = qs/2/kf;
[beta,gamma,Gm_inf] = getGm_ov(rs);
Gm = (beta*x.^2 + gamma*x.^4)./(1+gamma*x.^4/Gm_inf);
[alpha,Gp_inf] = getGp_ov(rs);
Gp = (1+alpha)*x.^2 ./(1 + (1+alpha)*x.^2/Gp_inf);
vf = (1-Gp).* vq;
Ixc= Gm .* vq;
%#**************************

qI = qs.^2 .* Ixc.^2;
qvq= qs.^2 .* vq.^2; 
qvf= qs.^2 .* vf.^2;


Sig_PP   = zeros(Nfreq,Nk);
ReSig_PP = zeros(Nfreq,Nk);
ReSig_RP = zeros(Nfreq,Nk);
ImSig_RP = zeros(Nfreq,Nk);
ReSig_SF = zeros(Nfreq,Nk);
ImSig_SF = zeros(Nfreq,Nk);
ReSig_LH = zeros(Nfreq,Nk);
ImSig_LH = zeros(Nfreq,Nk);
ReSig_OK = zeros(Nfreq,Nk);
ImSig_OK = zeros(Nfreq,Nk);
SigXs    = zeros(Nk,1);
SigXXs   = zeros(Nk,1);

for ik = 1:Nk;
  printf("doing ik=%f out of Nk=%f \n",ik,Nk);
  k = ks(ik);
  
  for iw = 1:Nfreq;
    w = ws(iw);
    
    SigPP_u  = zeros(Nu,1);
    ISigSF_u = zeros(Nu,1);
    ISigRP_u = zeros(Nu,1);
    ISigLH_u = zeros(Nu,1);
    ISigOK_u = zeros(Nu,1);
    SigX_u   = zeros(Nu,1);

    for iu = 1:Nu;
      u = us(iu);
      xi = (k^2+qs.^2-2*k*qs*u)/2 - ef;

      wshift = sign(xi).*(w-xi);
      x1 = qs/kf + ( wshift + I*eta)/ef./(qs/kf);
      x2 = qs/kf - ( wshift + I*eta)/ef./(qs/kf);
      f1 = (1 - x1.^2/4) .* log(  (x1+2)./(x1-2)  );
      f2 = (1 - x2.^2/4) .* log(  (x2+2)./(x2-2)  );
      u_qw = 0.5 + (f1 + f2)./(4*qs/kf); %# u -> 1 (q ->0, w=0)
      chi0 = chiP * u_qw; %# im(chi0) is anti-symmetric

      chiC = chi0 ./(1 -vq .* chi0);
      chiF = chi0 ./(1 -vf .* chi0);
      chiS = chi0 ./(1 -Ixc.* chi0);
      
      SigPP_u(iu)  = trapz(qs, 1./(wq .* (w-xi-wq.*sign(xi) + I*eta))); 
      ISigSF_u(iu) = trapz(qs, qI  .* imag(chiS) .* (wshift>0) );
      ISigRP_u(iu) = trapz(qs, qvq .* imag(chiC) .* (wshift>0) );
      ISigLH_u(iu) = trapz(qs, qvq .* imag(chiF) .* (wshift>0) );
      ISigOK_u(iu) = trapz(qs, qvf .* imag(chiF) .* (wshift>0) );
      SigX_u(iu)   = trapz(qs, (xi<0) );
    endfor; % u-loop
    
    ImSig_SF(iw,ik) += 3/(2*pi)^2 * trapz(us, ISigSF_u);
    ImSig_RP(iw,ik) += 1/(2*pi)^2 * trapz(us, ISigRP_u);
    ImSig_LH(iw,ik) += 1/(2*pi)^2 * trapz(us, ISigLH_u);
    ImSig_OK(iw,ik) += 1/(2*pi)^2 * trapz(us, ISigOK_u);
    Sig_PP(iw,ik)   += wp^2/(2*pi)* trapz(us, SigPP_u );

    if( iw == 1);
      SigXs(ik) += -1/pi * trapz(us,SigX_u); 
    endif;
  endfor; % w-loop

  %# bare exchange (Ashcroft Ch. 17, Eq. [17.19]):
  FX = 0.5 + (kf^2-k^2)/4/k/kf * log(abs( (kf+k)/(kf-k) ));
  if( abs(k-kf)  < 0.00001);
    FX = 0.5;
  elseif( abs(k) < 0.00001);
    FX = 1;
  endif;
  SigXXs(ik) = -2*kf/pi * FX;
  
  etaKK = 0.001/27.21;
  Sig_PP(:,ik)  += SigXXs(ik);
  ReSig_PP(:,ik) = getReSig(ws,-imag(Sig_PP(:,ik)),etaKK) + SigXXs(ik);
  ReSig_RP(:,ik) = getReSig(ws,-ImSig_RP(:,ik),etaKK) + SigXXs(ik);
  ReSig_LH(:,ik) = getReSig(ws,-ImSig_LH(:,ik),etaKK) + SigXXs(ik);
  ReSig_OK(:,ik) = getReSig(ws,-ImSig_OK(:,ik),etaKK) + SigXXs(ik);
  ReSig_SF(:,ik) = getReSig(ws,-ImSig_SF(:,ik),etaKK);


endfor; %# ik-loop

%# get dmu:
[a,ikf] = min( abs( ks-kf ));
printf("using k=%f * kf as Fermi wave vector \n",ks(ikf)/kf);
dmu_PP = spline(ws,real(Sig_PP(:,ikf)),0);
dmu_PP2= spline(ws,ReSig_PP(:,ikf),0);
dmu_RP = spline(ws,ReSig_RP(:,ikf),0);
dmu_LH = spline(ws,ReSig_LH(:,ikf),0);
dmu_OK = spline(ws,ReSig_OK(:,ikf),0);

ReSig_SF += ReSig_OK;
ImSig_SF += ImSig_OK;
dmu_SF    = spline(ws,ReSig_SF(:,ikf),0);

printf("rs = %f \n",rs);
printf("chemical potential mu = %f ry \n",dmu_PP*2);
printf("chemical potential mu2= %f ry \n",dmu_PP2*2);
printf("chemical potential rpa= %f ry \n",dmu_RP*2);
printf("chemical potential LH=  %f ry \n",dmu_LH*2);
printf("chemical potential OK=  %f ry \n",dmu_OK*2);
printf("chemical potential SF=  %f ry \n",dmu_SF*2);
printf(" \n");

DSig_PP = zeros(Nk,1);
DSig_SF = zeros(Nk,1);
DSig_OK = zeros(Nk,1);
DSig_LH = zeros(Nk,1);
DSig_RP = zeros(Nk,1);
for ik=1:Nk;
  ek = (ks(ik).^2/2 - ef);
  DSig_PP(ik) = interp1(ws,ReSig_PP(:,ik),ek);
  DSig_SF(ik) = interp1(ws,ReSig_SF(:,ik)-ReSig_OK(:,ik),ek);
  DSig_OK(ik) = interp1(ws,ReSig_OK(:,ik),ek);
  DSig_LH(ik) = interp1(ws,ReSig_LH(:,ik),ek);
  DSig_RP(ik) = interp1(ws,ReSig_RP(:,ik),ek);
endfor;

wss = [-5:0.01:1]'/27.21;
Nfreqs = length(wss);
A_PP = zeros(Nfreqs,Nk);
A_PP2= zeros(Nfreqs,Nk);
A_RP = zeros(Nfreqs,Nk);
A_SF = zeros(Nfreqs,Nk);
A_LH = zeros(Nfreqs,Nk);
A_OK = zeros(Nfreqs,Nk);
for ii = 1:Nk;
  xi = ks(ii)^2/2 -ef;
  ReSig = interp1(ws,real(Sig_PP(:,ii)),wss);
  ImSig = interp1(ws,imag(Sig_PP(:,ii)),wss);
  A_PP(:,ii) = abs(ImSig)./( (wss-xi-ReSig+dmu_PP).^2 + ImSig.^2);

  ReSig = interp1(ws,ReSig_PP(:,ii),wss);
  A_PP2(:,ii) = abs(ImSig)./( (wss-xi-ReSig+dmu_PP2).^2 + ImSig.^2);

  ReSig = interp1(ws,ReSig_RP(:,ii),wss);
  ImSig = interp1(ws,ImSig_RP(:,ii),wss);
  A_RP(:,ii) = abs(ImSig)./( (wss-xi-ReSig+dmu_RP).^2 + ImSig.^2);

  ReSig = interp1(ws,ReSig_SF(:,ii),wss);
  ImSig = interp1(ws,ImSig_SF(:,ii),wss);
  A_SF(:,ii) = abs(ImSig)./( (wss-xi-ReSig+dmu_SF).^2 + ImSig.^2);

  ReSig = interp1(ws,ReSig_LH(:,ii),wss);
  ImSig = interp1(ws,ImSig_LH(:,ii),wss);
  A_LH(:,ii) = abs(ImSig)./( (wss-xi-ReSig+dmu_LH).^2 + ImSig.^2);

  ReSig = interp1(ws,ReSig_OK(:,ii),wss);
  ImSig = interp1(ws,ImSig_OK(:,ii),wss);
  A_OK(:,ii) = abs(ImSig)./( (wss-xi-ReSig+dmu_OK).^2 + ImSig.^2);
endfor;
A_PP  /= pi;
A_PP2 /= pi;
A_RP  /= pi;
A_SF  /= pi;
A_LH  /= pi;
A_OK  /= pi;

D2Sig_PP = zeros(Nk,1);
D2Sig_SF = zeros(Nk,1);
D2Sig_OK = zeros(Nk,1);
D2Sig_LH = zeros(Nk,1);
D2Sig_RP = zeros(Nk,1);
for ik=1:Nk;
  [a,b] = max(A_PP(:,ik) );
  D2Sig_PP(ik) = interp1(ws,ReSig_PP(:,ik),wss(b));
  
  [a,b] = max(A_SF(:,ik) );
  D2Sig_SF(ik) = interp1(ws,ReSig_SF(:,ik)-ReSig_OK(:,ik),wss(b));
  
  [a,b] = max(A_OK(:,ik) );
  D2Sig_OK(ik) = interp1(ws,ReSig_OK(:,ik),wss(b));
  
  [a,b] = max(A_LH(:,ik) );
  D2Sig_LH(ik) = interp1(ws,ReSig_LH(:,ik),wss(b));
  
  [a,b] = max(A_RP(:,ik) );
  D2Sig_RP(ik) = interp1(ws,ReSig_RP(:,ik),wss(b));
endfor;


printf("qp-energies: \n");
[a,b] = max(A_PP(:,1));
printf("PP: %f ry \n",wss(b)*2);
[a,b] = max(A_PP2(:,1));
printf("PP2: %f ry \n",wss(b)*2);
[a,b] = max(A_RP(:,1));
printf("rpa: %f ry \n",wss(b)*2);
[a,b] = max(A_SF(:,1));
printf("SF: %f ry \n",wss(b)*2);
[a,b] = max(A_LH(:,1));
printf("LH: %f ry \n",wss(b)*2);
printf("H: %f \n",-ef*2);

Sig_LH = ReSig_LH + I*ImSig_LH;
Sig_RP = ReSig_RP + I*ImSig_RP;
Sig_OK = ReSig_OK + I*ImSig_OK;
Sig_SF = ReSig_SF + I*ImSig_SF;

save output ws Sig_PP Sig_RP Sig_OK Sig_LH Sig_SF rs ks

more on;