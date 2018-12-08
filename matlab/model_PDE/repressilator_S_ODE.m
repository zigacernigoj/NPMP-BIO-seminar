function  [dmA, dmB, dmC, dA, dB, dC, dS_i, dS_e] = repressilator_S_ODE(CELLS, mA, mB, mC, A, B, C, S_i, S_e, alpha, alpha0, Kd, beta, delta_m, delta_p, n, kS0, kS1, kSe, kappa, eta)
    dmA = CELLS .*(alpha./(1 + (C./Kd).^n) + alpha0 - delta_m * mA);
    dmB = CELLS .*(alpha./(1 + (A./Kd).^n) + alpha0 - delta_m * mB);
    dmC = CELLS .*(alpha./(1 + (B./Kd).^n) + alpha0 - delta_m * mC + (kappa * S_i)./(1 + S_i));
    
    dA = CELLS .*(beta * mA - delta_p * A);
    dB = CELLS .*(beta * mB - delta_p * B);
    dC = CELLS .*(beta * mC - delta_p * C);
    
    dS_i = CELLS .* (- kS0 * S_i + kS1 * A - eta * (S_i - S_e));
    dS_e = - kSe * S_e + CELLS .*(eta * (S_i - S_e));
end

    