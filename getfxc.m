function [fxc] = getfxc(n);

alpha = 3/4 * (3/2/pi)^(2/3);
a =  0.0311;
b = -0.0480;
c =  0.0020;
d = -0.0116;
A = -0.1423;
B =  1.0529;
C =  0.3334;

rs = (3/4/pi)^(1/3) * n.^(-1/3);

%# get first derivate: excp = d exc / dn = (d exc / drs) * drs/dn
rsp = (3/4/pi)^(1/3) * (-1/3.) * n.^(-4/3);
%# calculate dexc/drs:
dedrs = alpha./rs.^2;
if(rs < 1.0)
  dedrs += a./rs + c*(log(rs)+1)+ d;
else;
  dedrs += -A./(1+B*sqrt(rs)+C*rs).^2 .*( B/2./sqrt(rs)+C );
endif;

excp = dedrs .* rsp;

%# second derivative: excpp = d^2 exc / dn^2 = (d^2 exc / drs^2) * (drs/dn)^2 + dexc/drs * (d^2 rs/dn^2)
rspp = (3/4/pi)^(1/3) * (4/9.) .* n.^(-7/3);
%# get second derivative: exc
ddedrs = -2*alpha./rs.^3;
if(rs<1);
  ddedrs += -a./rs.^2 + C./rs;
else;
  ddedrs += ( 2*(B/2./sqrt(rs) +C ).^2 + B/4./rs.^(3/2) .* (1+B*sqrt(rs) + C*rs) ) *A./(1+B*sqrt(rs)+C*rs).^3;
endif;
excpp = ddedrs .* rsp.^2 + dedrs .* rspp;

fxc = 2*excp + n.*excpp;
endfunction;