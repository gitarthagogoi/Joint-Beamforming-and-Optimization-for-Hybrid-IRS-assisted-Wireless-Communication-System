function [W, Theta, Rsum] = hybrid_RIS_precoding_zetaxi(M,K,N,Ps_max,Pr_cap,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G,Ka_input)
% HYBRID_RIS_PRECODING_ZETAXI
% Hybrid RIS using zeta/xi selection (paper-inspired) + AO/WMMSE loop.
%
% Usage:
% [W, Theta, Rsum] = hybrid_RIS_precoding_zetaxi(M,K,N,Ps_max,Pr_cap,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G,Ka_input)
%
% Inputs:
%  - M,K,N : system dims (BS antennas, users, RIS elements)
%  - Ps_max, Pr_cap : BS and RIS amplifier power (Watts)
%  - sigma2, sigmar2 : noise powers
%  - eta_k : Kx1 user weights
%  - Theta : initial RIS diag (NxN or Nx1 diag-vector)
%  - W     : initial precoder (MxK matrix or vectorized length M*K)
%  - h_k,f_k,G : channels (as in your code): f_k is KxN, G is N x M
%  - Ka_input : number of active elements to pick (optional)
%
% Outputs:
%  - W : final BS precoder (MxK)
%  - Theta : final hybrid RIS diagonal matrix (NxN)
%  - Rsum : achieved sum-rate (scalar)
%
% Notes:
%  - This selects active set once from the passive warm-start (you can adapt to reselect each AO iter)
%  - Implemented with per-element KxK matrices (users dimension) to approximate An,Bn,Cn
%  - Regularized solves and eig fallbacks included.

    if nargin < 13
        Ka = max(1, round(0.05 * N));
    else
        Ka = Ka_input;
        Ka = max(1, min(N, round(Ka)));
    end

    eps_small = 1e-12;
    maxIter = 30;
    Rsum_hist = zeros(2*maxIter,1);

    % ----------------------
    % 1) Passive warm-start (get good phases)
    % ----------------------
    [W_pass, Theta_pass, Rsum_pass] = passive_RIS_precoding(M,K,N,Ps_max,sigma2,eta_k,Theta,W,h_k,f_k,G);
    W = W_pass;
    if isvector(Theta_pass), theta_vec = Theta_pass(:); else theta_vec = diag(Theta_pass); end
    Theta = diag(theta_vec);

    % Convert W to MxK matrix if needed
    if ismatrix(W) && all(size(W) == [M, K])
        W_mat = W;
    elseif isvector(W) && numel(W) == M*K
        W_mat = reshape(W, M, K);
    else
        try
            wtmp = w_k_generate(K,M,W);
            W_mat = w_k2W(K,M,wtmp);
        catch
            error('hybrid_RIS_precoding_zetaxi: cannot convert W to MxK. Provide correct shape or helper functions.');
        end
    end

    % ----------------------
    % 2) Compute zeta_n and xi_n for each RIS element
    %    We follow the spirit of the paper using KxK matrices:
    %      - rn := f_k(:,n)  (K x 1)
    %      - tn := (G(n,:) * W_mat).' (K x 1)  (proxy for BS->element->users)
    %    Build sumAi = sum_{i != n} alpha_i * (rn_i * tn_i') -> K x K
    %    An = I + rho * sumAi * sumAi'
    %    Bn = rho * rn * ( (tn' * tn) * rn' )
    %    Cn = rho * rn * ( tn' * sumAi' )
    %    gamma_n = trace(An \ Bn)
    %    lambda_n = principal eigenvalue of En^{-1} Cn with En = An * (I + An\Bn)
    % ----------------------

    rho = Ps_max / max(sigma2, eps_small);
    zeta = zeros(N,1);
    xi   = zeros(N,1);

    % Precompute per-element tn and rn
    t_mat = zeros(K,N); % each column n is tn (Kx1) flattened
    for n = 1:N
        % G(n,:) 1 x M, W_mat M x K -> 1 x K, transpose -> Kx1
        t_mat(:,n) = (G(n,:) * W_mat).'; 
    end
    r_mat = f_k; % K x N

    % Precompute outer products maybe heavy; we compute per n loop
    for n = 1:N
        % build sumAi (K x K)
        sumAi = zeros(K,K);
        for i = 1:N
            if i == n, continue; end
            ai = theta_vec(i);             % complex coefficient (phase * amp)
            ri = r_mat(:, i);              % K x 1
            ti = t_mat(:, i);              % K x 1
            sumAi = sumAi + ai * (ri * ti'); % K x K
        end

        An = eye(K) + rho * (sumAi * sumAi');    % K x K
        rn = r_mat(:, n);                         % K x 1
        tn = t_mat(:, n);                         % K x 1

        % Bn, Cn
        % tn' * tn is scalar
        Bn = rho * ( rn * ( (tn' * tn) * rn' ) );   % K x K (rank-1)
        Cn = rho * ( rn * ( tn' * (sumAi') ) );     % K x K

        % regularize An
        An_reg = An + eps_small * eye(K);

        % gamma_n approx as trace(An^{-1} Bn)
        gamma_n = real(trace(An_reg \ (Bn + eps_small * eye(K))));

        % En and principal eigenvalue
        Dn = eye(K) + (An_reg \ Bn);
        En = An_reg * Dn + eps_small * eye(K);

        % compute lamMat = En^{-1} * Cn and principal eigenvalue
        try
            lamMat = En \ (Cn + eps_small*eye(K));
            ev = eig(lamMat);
            if isempty(ev)
                lambda_n = 0;
            else
                [~, idxMax] = max(abs(ev));
                lambda_n = ev(idxMax);
            end
        catch
            % fallback: scalar proxy using trace
            lamMat = En \ (Cn + eps_small*eye(K));
            lambda_n = trace(lamMat);
        end

        zeta(n) = abs(lambda_n) * sqrt(max(0, gamma_n));
        xi(n)   = sigma2 + Ps_max * (norm(tn)^2);
    end

    % select top-Ka by zeta/xi
    scores = zeta ./ (xi + 1e-16);
    [~, idx_sorted] = sort(scores, 'descend');
    active_idx = idx_sorted(1:min(Ka, length(idx_sorted)));
    passive_idx = setdiff(1:N, active_idx);

    % ----------------------
    % 3) AO / WMMSE loop (unchanged) — use cvx_solve_W and cvx_solve_theta as before
    %    enforce passive elements amplitude = 1 after theta update
    % ----------------------
    % Ensure Theta is diag matrix
    Theta = diag(theta_vec);

    for Q = 1:maxIter
        % update W
        w_k = w_k_generate(K,M,W);  % helper in your code
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

        % update Theta
        eps_k = eps_update(K,M,N,Rho_k,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);
        [nu, Lam] = nu_Lam_generate(K,M,N,Rho_k,eta_k,eps_k,h_k,f_k,G,w_k,sigmar2);

        theta_init = diag(Theta);
        try
            theta_new = cvx_solve_theta(N,K,M,theta_init,nu,Lam,w_k,G,Pr_cap,sigmar2);
        catch ME
            rethrow(ME);
        end

        % enforce passive amplitude = 1
        amp = abs(theta_new);
        phase = angle(theta_new);
        amp(passive_idx) = 1;
        theta_new = amp .* exp(1j * phase);

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

    % ----------------------
    % 4) Return best sum-rate
    % ----------------------
    Rsum = max(Rsum_hist(1:2*Q));
    % outputs W and Theta are already updated
end
