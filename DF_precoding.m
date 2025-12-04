function [Rsum_DF] = DF_precoding(K, Ps_max, Pr_max, sigma2, f_k, G)
% Simple Decode-and-Forward Rate Computation
% Inputs:
%   K       - number of users
%   Ps_max  - BS transmit power
%   Pr_max  - RIS (relay) transmit power
%   sigma2  - noise power
%   f_k     - RIS->User channels (KxN)
%   G       - BS->RIS channel (NxM)
%
% Output:
%   Rsum_DF - total DF sum rate for K users (bps/Hz)

Rsum_DF = 0;
gain_BR = mean(abs(G(:)).^2);  % BS->RIS
for k = 1:K
    gain_RU = mean(abs(f_k(k,:)).^2);  % RIS->User
    gamma1 = Ps_max * gain_BR / sigma2;
    gamma2 = Pr_max * gain_RU / sigma2;
    Rsum_DF = Rsum_DF + 0.5 * log2(1 + min(gamma1, gamma2));
end
end
