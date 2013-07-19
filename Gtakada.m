function [Ax,Ac,Bx,Bc,kappaF,g0] = Gtakada(rs);

  mmax = 30;
  alpha  = (4/9/pi)^(1/3);
  lambda = sqrt(4*alpha*rs/pi);
  a1 = 12.05;
  a2 = 4.254;
  a3 = 1.363;

  kf = (9*pi/4)^(1/3)/rs; %# see Ashcroft (a.u.)
  Nf = kf/pi^2;

  kappaF = 1 - lambda^2/4 * (1 + 0.07671*lambda^2*( (1+a1*lambda)^2 + 4/3*a2*lambda^2*(1+7/8*a1*lambda)+3/2*a3*lambda^3*(1+8/9*a1*lambda)) / (1+a1*lambda+a2*lambda^2+a3*lambda^3)^2 );

  chiF   = 1 - lambda^2/4 * (1 + lambda^2/8 * (log(lambda^2/(lambda^2+0.990)) - (1.122+1.222*lambda^2)/(1+0.533*lambda^2+0.184*lambda^4)));

  g0 = 0;
  for m = 0:mmax;
    g0 += lambda^(2*m)/factorial(m)/factorial(m+1);
  endfor;
  g0 = 1/g0^2;

  Ax = -kf^2/8/pi/Nf*(kappaF + chiF -2);
  Ac = -kf^2/8/pi/Nf*(kappaF - chiF );
  Bx = 6*Ax/(1+g0);
  Bc = 2*Ac/(1-g0);

#  printf("test: \n");
#  printf("q -> inf: G(+)= %f and %f \n",Ax/Bx+Ac/Bc,2/3-g0/3);
#  printf("q -> inf: G(-)= %f and %f \n",Ax/Bx-Ac/Bc,2/3*g0-1/3);
#  printf("q ->   0: G(+)= %f and %f \n",1+4*pi*Nf*(Ax+Ac)/kf^2,kappaF);
#  printf("q ->   0: G(-)= %f and %f \n",1+4*pi*Nf*(Ax-Ac)/kf^2,chiF);
  
endfunction;