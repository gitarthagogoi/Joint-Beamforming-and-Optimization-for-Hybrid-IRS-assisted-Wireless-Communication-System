function [W, Theta, Rsum] = hybrid_RIS_precoding_simple(M,K,N,Ps_max,Pr_cap,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G)
% HYBRID_RIS_PRECODING
%   Hybrid RIS precoding based on passive optimization + amplifier allocation.
%   Pr_cap is total RIS amplifier budget in Watts (e.g. 0.01 for 10 dBm).
%
% Inputs/Outputs match your active/passive functions.

iteration = 30;

% --- Step 1: run passive optimization first to get good phases ---
[W_passive, Theta_passive, Rsum_passive] = passive_RIS_precoding(M,K,N,Ps_max,sigma2,eta_k,Theta,W,h_k,f_k,G);

% extract passive theta (N×1)
theta = diag(Theta_passive);

% --- Step 2: choose Ka active elements (e.g. strongest ones) ---
Ka = max(1, round(0.05*N));  % 5% active
strength = sqrt(sum(abs(G).^2,2)) .* sum(abs(f_k),1).';  % cascaded strength metric
[~, idx] = sort(strength,'descend');
active_idx = idx(1:Ka);

% --- Step 3: allocate amplifier power ---
alpha = ones(N,1);
P_per = Pr_cap / Ka;

% estimate input power to each active element
avg_gain = mean(abs(G).^2,2);   % N×1
input_power = avg_gain * (Ps_max/M) + sigma2;
denom_vec = input_power(active_idx);

alpha(active_idx) = sqrt(1 + P_per ./ (denom_vec + eps));

% scale if exceeding total Pr_cap
P_ris_actual = sum((alpha(active_idx).^2 - 1) .* denom_vec);
if P_ris_actual > Pr_cap
    scale = sqrt(Pr_cap / (P_ris_actual + eps));
    alpha(active_idx) = 1 + (alpha(active_idx)-1)*scale;
end

% --- Step 4: build hybrid Theta ---
Theta = diag(alpha .* theta);

% --- Step 5: recompute effective SINR and sum-rate ---
w_k = w_k_generate(K,M,W_passive);
[Rsum,~] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

end
