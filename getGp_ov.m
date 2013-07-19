function [alpha,Gp_inf] = getGp_ov(rs);

  alpha  = 0.0335/2 +  0.02*rs/3/(0.1+rs)^2 * (1+rs/(0.1+rs) );
  alpha *= (2*pi/3)^(2/3)*rs;
  g0 = 32/(8+3*rs)^2;
  Gp_inf = 2/3 * ( 1- g0 );

endfunction;