function [W, Theta, Rsum] = hybrid_RIS_precoding_phasealignment(M,K,N,Ps_max,Pr_cap,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G,Ka_input)
% HYBRID_RIS_PRECODING_PHASEALIGNMENT
% Hybrid RIS with phase-alignment based active selection (one-shot selection,
% then AO/WMMSE loop unchanged).
%
% Inputs/Outputs: same conventions as your other hybrid functions.
% Ka_input: optional number of active elements (default 5% of N).

    if nargin < 13
        Ka = max(1, round(0.05 * N));
    else
        Ka = Ka_input;
        Ka = max(1, min(N, round(Ka)));
    end

    maxIter = 30;
    Rsum_hist = zeros(2*maxIter,1);

    % -- ensure Theta is diagonal vector form
    if isvector(Theta)
        theta_vec = Theta(:);
    else
        theta_vec = diag(Theta);
    end
    Theta = diag(theta_vec);

    % -------------------------
    % Build W_mat (M x K) robustly
    % -------------------------
    % Accept:
    %  - W as MxK matrix
    %  - W as vector length M*K (column-major)
    %  - fallback to helper conversions if available
    if ismatrix(W) && all(size(W) == [M, K])
        W_mat = W;
    elseif isvector(W) && numel(W) == M*K
        W_mat = reshape(W, M, K);
    else
        % try to use your helper conversion (if present)
        try
            w_k_tmp = w_k_generate(K,M,W);     % may error depending on W format
            W_mat = w_k2W(K,M,w_k_tmp);        % returns MxK
        catch
            error('Cannot convert input W to MxK. Provide W as MxK matrix or vector length M*K.');
        end
    end

    % --- Phase-alignment selection (one-shot, using current Theta and W_mat)
    % Compute per-element score: abs( f_k(:,n)' * (theta_n * (G(n,:) * W_mat).' ) )
    score = zeros(N,1);
    theta_diag = diag(Theta); % N x 1
    for n = 1:N
        % G(n,:) is 1 x M; multiply by W_mat (M x K) -> 1 x K row vector of complex contributions
        tx_row = G(n,:) * W_mat;    % 1 x K
        tx_contrib = tx_row.';      % K x 1 (per-user complex contributions)
        fk_n = f_k(:, n);           % K x 1
        % scalar projection onto user channels after applying RIS phase theta_n
        proj = fk_n' * (theta_diag(n) * tx_contrib); % 1x1 complex
        score(n) = abs(proj) + 0;   % keep real nonnegative scalar
    end

    [~, idx_sorted] = sort(score, 'descend');
    active_idx = idx_sorted(1:min(Ka,length(idx_sorted)));
    passive_idx = setdiff(1:N, active_idx);

    % --- Ensure Theta is diagonal matrix for AO loop
    Theta = diag(theta_diag);

    % --- AO / WMMSE loop (identical structure to other hybrid functions)
    for Q = 1:maxIter
        % Update W (BS precoder)
        w_k = w_k_generate(K,M,W);  % your helper (may accept vectorized W)
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

        % Update Theta via your CVX solver
        eps_k = eps_update(K,M,N,Rho_k,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);
        [nu, Lam] = nu_Lam_generate(K,M,N,Rho_k,eta_k,eps_k,h_k,f_k,G,w_k,sigmar2);

        theta_init = diag(Theta);
        try
            theta_new = cvx_solve_theta(N,K,M,theta_init,nu,Lam,w_k,G,Pr_cap,sigmar2);
        catch ME
            rethrow(ME);
        end

        % enforce passive elements amplitude = 1, allow amplification only on active_idx
        amp = abs(theta_new);
        ph  = angle(theta_new);
        amp(passive_idx) = 1;
        theta_new = amp .* exp(1j*ph);

        Theta = diag(theta_new);

        [Rsum_hist(2*Q), gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

        % convergence check
        if Q > 1
            prev = Rsum_hist(2*Q-2);
            curr = Rsum_hist(2*Q);
            if prev > 0 && (curr - prev)/prev < 0.001
                break;
            end
        end
    end

    % Return best Rsum over history
    Rsum = max(Rsum_hist(1:2*Q));
    % ensure outputs in expected shapes
    % W is already updated by AO loop (MxK)
end
