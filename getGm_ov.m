function [beta,gamma,Gm_inf] = getGm_ov(rs);

  chi_o_chip = 1.17 + 0.029*(rs-1)^2 + 0.175*log(rs);
  beta = 1/0.166/rs * (1-1/chi_o_chip);

  g0 = 32/(8+3*rs)^2;
  Gm_inf = 1/3 * ( 4*g0 - 1);
  gamma = beta*Gm_inf/(-2*Gm_inf + beta);

endfunction;