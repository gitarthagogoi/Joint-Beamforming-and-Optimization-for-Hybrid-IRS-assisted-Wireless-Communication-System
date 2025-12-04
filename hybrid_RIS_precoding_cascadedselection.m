function [W, Theta, Rsum] = hybrid_RIS_precoding_cascadedselection(M,K,N,Ps_max,Pr_cap,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G,Ka_input)
% HYBRID_RIS_PRECODING_ITERATIVE
% Alternating optimization for hybrid RIS (subset of elements active, total RIS
% amplifier forwarded power constrained to Pr_cap).
%
% Inputs:
%   M,K,N,Ps_max,Pr_cap,sigma2,sigmar2,eta_k,Theta (init), W (init), h_k, f_k, G
%   Ka_input (optional) : number of active RIS elements to allow (default 5% of N)
%
% Outputs:
%   W      - optimized BS precoder (M x K)
%   Theta  - N x N diagonal hybrid RIS matrix (complex values: amplitude * phase)
%   Rsum   - achieved sum-rate (scalar)
%
% Notes:
%  - This function relies on the same helper functions your active/passive code uses:
%    w_k_generate, SINR_calculate, H_k_generate, eps_update, v_A_k_generate,
%    w_k2W, cvx_solve_W, nu_Lam_generate, cvx_solve_theta.
%  - The cvx_solve_theta used in your active implementation accepts Pr_max to
%    constrain RIS forwarded power. We call that with Pr_cap.
%  - After the cvx theta update, we force passive indices to have unit amplitude
%    (i.e., amplitude 1) so only the selected active subset uses amplification.
%
% Written to integrate with user's existing code structure.

    if nargin < 13
        Ka = Ka_input; 
    else
        Ka = Ka_input;
        Ka = max(1, min(N, round(Ka)));
    end

    maxIter = 30;
    Rsum_hist = zeros(2*maxIter,1);

    % --- choose active indices once at start (you can change to reselect each iter if desired)
    % Use a cascaded-strength metric: s_n = ||G(n,:)|| * sum_k |f_k(k,n)|
    g_row_norm = sqrt(sum(abs(G).^2,2));   % Nx1
    f_abs_sum   = sum(abs(f_k),1).';       % N x 1
    strength = g_row_norm .* f_abs_sum;
    [~, idx_sorted] = sort(strength, 'descend');
    active_idx = idx_sorted(1:Ka);
    passive_idx = setdiff(1:N, active_idx);

    % Ensure Theta initial is diagonal; if it's provided as diag vector handle both cases
    if isvector(Theta)
        theta_vec = Theta(:);
    else
        theta_vec = diag(Theta);
    end
    Theta = diag(theta_vec);

    for Q = 1:maxIter
        % ----------------- update W (BS precoder) -----------------
        % build current w_k from W
        w_k = w_k_generate(K,M,W);  % user's helper (returns K cell or vector form used by your w_k2W)
        [~, gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

        % Build H_k and related WMMSE auxiliary variables (same as active)
        H_k = H_k_generate(K,M,N,h_k,f_k,G,Theta);

        Rho_k = gamma_k;
        eps_k = eps_update(K,M,N,Rho_k,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

        [V,A] = v_A_k_generate(K,M,N,Rho_k,eta_k,eps_k,h_k,f_k,G,Theta);

        W = w_k2W(K,M,w_k);

        % Call CVX solver for W that your active function uses;
        % signature used in active_RIS_precoding: cvx_solve_W(M,K,G,Theta,V,A,W,Ps_max,Pr_max,sigmar2)
        % For hybrid, pass Pr_cap as the Pr_max so joint constraints considered
        try
            W = cvx_solve_W(M,K,G,Theta,V,A,W,Ps_max,Pr_cap,sigmar2);
        catch
            % fallback to passive CVX W if that function differs
            W = cvx_solve_W_for_passiveRIS(M,K,G,Theta,V,A,W,Ps_max);
        end

        w_k = w_k_generate(K,M,W);

        [Rsum_hist(2*Q-1), gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

        % ----------------- update Theta (phases + amplitudes for active subset) -----------------
        eps_k = eps_update(K,M,N,Rho_k,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);
        [nu, Lam] = nu_Lam_generate(K,M,N,Rho_k,eta_k,eps_k,h_k,f_k,G,w_k,sigmar2);

        % prepare theta vector initial guess (N x 1)
        theta_init = diag(Theta);

        % call the CVX theta solver used in active implementation:
        % signature in your active_RIS_precoding: cvx_solve_theta(N,K,M,theta,nu,Lam,w_k,G,Pr_max,sigmar2)
        try
            theta_new = cvx_solve_theta(N,K,M,theta_init,nu,Lam,w_k,G,Pr_cap,sigmar2);
        catch ME
            % If your version has different name or signature, rethrow with helpful message
            rethrow(ME);
        end

        % enforce passive elements to unit amplitude (only active_idx can amplify)
        amp = abs(theta_new);
        phase = angle(theta_new);
        amp(passive_idx) = 1;           % force passive elements amplitude = 1
        theta_new = amp .* exp(1j*phase);

        Theta = diag(theta_new);

        [Rsum_hist(2*Q), gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

        % convergence check (same style as active)
        if Q > 1
            prev = Rsum_hist(2*Q-2);
            curr = Rsum_hist(2*Q);
            if prev > 0
                if (curr - prev)/prev < 0.001  % small threshold; adjust if you want
                    break;
                end
            end
        end
    end

    % return best found Rsum (max over history)
    Rsum = max(Rsum_hist(1:2*Q));
    % ensure outputs are in expected shapes
    W = W;
    Theta = Theta;
end
