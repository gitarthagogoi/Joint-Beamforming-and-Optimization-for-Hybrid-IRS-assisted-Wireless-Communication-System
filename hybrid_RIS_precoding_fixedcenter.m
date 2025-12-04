function [W, Theta, Rsum] = hybrid_RIS_precoding_fixedcenter(M,K,N,Ps_max,Pr_cap,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G,Ka_input)
% HYBRID_RIS_PRECODING_FIXEDCENTER
% Alternating optimization for hybrid RIS with fixed-center active elements.
%
% Inputs:
%   M,K,N,Ps_max,Pr_cap,sigma2,sigmar2,eta_k,Theta(init),W(init),h_k,f_k,G
%   Ka_input (optional) : number of active RIS elements (default 5% of N)
%
% Outputs:
%   W      - optimized BS precoder (M x K)
%   Theta  - N x N diagonal RIS matrix (hybrid: Ka active + others passive)
%   Rsum   - achieved sum-rate (scalar)
%
% Notes:
%  - Same AO/WMMSE logic as hybrid_RIS_precoding_cascadedselection.
%  - Difference: active_idx is chosen as a fixed block centered at N/2.

    if nargin < 13
        Ka = max(1, round(0.05 * N)); % default 5% of N
    else
        Ka = Ka_input;
        Ka = max(1, min(N, round(Ka)));
    end

    maxIter = 30;
    Rsum_hist = zeros(2*maxIter,1);

    % --- Fixed-center active index selection ---
    center = ceil(N/2);
    half = floor(Ka/2);
    startIdx = max(1, center - half);
    endIdx   = min(N, startIdx + Ka - 1);
    active_idx = startIdx:endIdx;
    passive_idx = setdiff(1:N, active_idx);

    % --- Ensure Theta is diagonal vector ---
    if isvector(Theta)
        theta_vec = Theta(:);
    else
        theta_vec = diag(Theta);
    end
    Theta = diag(theta_vec);

    for Q = 1:maxIter
        % ----------------- update W (BS precoder) -----------------
        w_k = w_k_generate(K,M,W);
        [~, gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

        H_k = H_k_generate(K,M,N,h_k,f_k,G,Theta);

        Rho_k = gamma_k;
        eps_k = eps_update(K,M,N,Rho_k,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

        [V,A] = v_A_k_generate(K,M,N,Rho_k,eta_k,eps_k,h_k,f_k,G,Theta);
        W = w_k2W(K,M,w_k);

        try
            W = cvx_solve_W(M,K,G,Theta,V,A,W,Ps_max,Pr_cap,sigmar2);
        catch
            W = cvx_solve_W_for_passiveRIS(M,K,G,Theta,V,A,W,Ps_max);
        end

        w_k = w_k_generate(K,M,W);
        [Rsum_hist(2*Q-1), gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

        % ----------------- update Theta -----------------
        eps_k = eps_update(K,M,N,Rho_k,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);
        [nu,Lam] = nu_Lam_generate(K,M,N,Rho_k,eta_k,eps_k,h_k,f_k,G,w_k,sigmar2);

        theta_init = diag(Theta);
        try
            theta_new = cvx_solve_theta(N,K,M,theta_init,nu,Lam,w_k,G,Pr_cap,sigmar2);
        catch ME
            rethrow(ME);
        end

        % enforce passive amplitude=1, active can amplify
        amp = abs(theta_new);
        phase = angle(theta_new);
        amp(passive_idx) = 1;
        theta_new = amp .* exp(1j*phase);

        Theta = diag(theta_new);

        [Rsum_hist(2*Q), gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

        % convergence check
        if Q > 1
            prev = Rsum_hist(2*Q-2);
            curr = Rsum_hist(2*Q);
            if prev > 0
                if (curr - prev)/prev < 0.001
                    break;
                end
            end
        end
    end

    Rsum = max(Rsum_hist(1:2*Q));
    W = W;
    Theta = Theta;
end
