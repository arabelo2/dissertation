function [E_n, P_n, Pi_n, phii_n] = Normalization(E, P, Pi, phii, c, delta_t, rho, Uo, K);

factor_phii = 1/c;
factor_Pi = factor_phii*delta_t/rho;
factor_P = factor_Pi/Uo;
%factor_E = factor_P*factor_P*(2*rho*c/ka);
%factor_E = (2/(K*rho))*delta_t^2/(Uo*c);
factor_E = 2*delta_t^2/(Uo*K*rho*c);

phii_n = factor_phii.*phii;
Pi_n = factor_Pi.*Pi;
P_n = factor_P.*P;
E_n = factor_E.*E;

% It does not seem it works as expected.
%phii_n = phii./max(abs(phii));
%Pi_n = Pi./max(abs(Pi));
%P_n = P./max(abs(P));
%E_n = E./max(abs(E));