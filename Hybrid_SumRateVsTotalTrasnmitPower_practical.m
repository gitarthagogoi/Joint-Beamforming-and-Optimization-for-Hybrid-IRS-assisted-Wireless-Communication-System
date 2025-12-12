tic

%% Parameters
user_X = 300;                 % distance setting (BS–RIS–Users geometry)
ExNumber = 2;                % Monte Carlo trials per power point (increase from 1)
M = 4;                        % BS transmit antennas
K = 4;                        % number of users
N = 512;                      % RIS elements

Ps_max_dBm = -10:5:40;        % sweep BS transmit power (dBm)
Ps_W = 10.^((Ps_max_dBm-30)/10);   % convert dBm to Watts

% Choose RIS amplifier budget setting
use_fixed_Pr = false;   % toggle between fixed or fraction
if use_fixed_Pr
    Pr_cap_fixed = 0.01;   % 10 dBm = 0.01 W fixed amplifier budget
end
frac_Pr = 0.01;             % fraction of BS power to give RIS if not using fixed (1% default)

Ka_input = 128;              % number of active elements for iterative hybrid

% Noise
sigma2 = 1e-10;              % check units (Watts)
sigmar2 = sigma2;            % amplifier noise (same order)

% Pathloss params
f_c = 5;                     % GHz
large_fading_AI = 2.2;
large_fading_DI = 2.2;

% weights
eta_k = ones(K,1);

% RIS practical phase quantization settings
Q = 4;                      % number of quantization bits (change as needed)
levels = 2^Q;               % number of discrete phase levels
phase_levels = (0:(levels-1)) * (2*pi/levels);  % allowed phase values (radians)

%% Preallocate results
nP = length(Ps_W);
Rsum_active   = zeros(nP, ExNumber);
Rsum_passive  = zeros(nP, ExNumber);
Rsum_random   = zeros(nP, ExNumber);
Rsum_noRIS    = zeros(nP, ExNumber);
Rsum_hybrid   = zeros(nP, ExNumber);
Rsum_hybrid_iter = zeros(nP, ExNumber);

%% Main sweep
% outer loop not parfor (parfor inside or parfor over b)
for a = 1:nP
    Ps = Ps_W(a);

    % choose RIS amplifier budget
    if use_fixed_Pr
        Pr_cap = Pr_cap_fixed;
    else
        Pr_cap = frac_Pr * Ps;   % use fraction of BS power (consistent)
    end

    fprintf('>>> Power Point %d/%d: Ps=%.1f dBm (%.6f W), Pr_cap=%.6f W\n',...
        a,nP,Ps_max_dBm(a),Ps,Pr_cap);

    for b = 1:ExNumber
        % reproducible RNG per trial:
        rng(1000 + a*1000 + b,'twister');

        % --- Generate channels ---
        [Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser] = Position_generate_2(K,user_X);
        [h_k,f_k,G] = Channel_generate2(K,N,M,large_fading_AI,large_fading_DI,...
                                       Dis_BStoRIS,Dis_BStoUser,Dis_RIStoUser,f_c);

        % --- Initializations ---

        % Create continuous RIS phases then quantize them to Q bits (mid-rise style)
        theta_cont = 2*pi*rand(N,1);                % continuous phases in [0,2pi)
        % Map to indices nearest quantization level (0,...,levels-1)
        indices = round(theta_cont * levels / (2*pi));  % nearest index
        indices = mod(indices, levels);                  % wrap so index==levels -> 0
        quant_phases = (2*pi/levels) * indices;         % quantized phases in [0,2pi)
        Theta0 = diag(exp(1j * quant_phases));         % quantized RIS diagonal (practical)

        % small sanity check (optional) - ensure quantized phases are among allowed set
        % chk = unique(mod(angle(diag(Theta0)),2*pi));
        % assert(all(ismember(round(chk,12), round(phase_levels,12))), 'Quantization mismatch');

        % --- Initialize W0 as M x K matrix and normalize to total power Ps ---
        W0_mat = exp(1j*2*pi*rand(M,K));         % random BS precoder initial phases
        % scale each column to unit norm then normalize whole matrix to power Ps
        for kk = 1:K
            W0_mat(:,kk) = W0_mat(:,kk)/norm(W0_mat(:,kk) + eps);
        end
        % scale so that sum of column norms squared = Ps
        current_power = sum(sum(abs(W0_mat).^2));
        if current_power == 0
            current_power = eps;
        end
        W0_mat = W0_mat * sqrt(Ps/current_power);
        % If your downstream functions expect vectorized W, reshape there or adapt accordingly.

        % --- No RIS baseline ---
        [~, Rsum_noRIS(a,b)] = NoRIS_precoding(M,K,N,Ps,sigma2,eta_k,W0_mat,h_k,f_k,G);

        % --- Random RIS (but use quantized Theta0 for consistency) ---
        [~, ~, Rsum_random(a,b)] = random_RIS_precoding(M,K,N,Ps,sigma2,eta_k,Theta0,W0_mat,h_k,f_k,G);

        % --- Passive RIS (use quantized Theta0 as initial/practical phases) ---
        [W_pass, Theta_pass, Rsum_passive(a,b)] = passive_RIS_precoding(M,K,N,Ps,sigma2,eta_k,Theta0,W0_mat,h_k,f_k,G);

        % Ensure Theta_pass stays quantized (some implementations return continuous Theta)
        % If passive_RIS_precoding returns continuous Theta, re-quantize it to hardware levels:
        if exist('Theta_pass','var') && ~isempty(Theta_pass)
            % extract diag angles, quantize again
            tp = mod(angle(diag(Theta_pass)), 2*pi);
            idx_tp = round(tp * levels / (2*pi));
            idx_tp = mod(idx_tp, levels);
            Theta_pass = diag(exp(1j*(2*pi/levels)*idx_tp));
        end

        % --- Active RIS (use Pr_cap) ---
        [~, ~, Rsum_active(a,b)] = active_RIS_precoding(M,K,N,Ps-Pr_cap, Pr_cap, sigma2,sigmar2,eta_k,Theta_pass,W_pass,h_k,f_k,G);

        % --- Iterative Hybrid RIS (AO, Ka active) ---
        [~, ~, Rsum_hybrid_iter(a,b)] = hybrid_RIS_precoding_cascadedselection(M,K,N,Ps,Pr_cap,sigma2,sigmar2,eta_k,Theta_pass,W_pass,h_k,f_k,G,Ka_input);

    end
end

%% Average across trials
Rsum_active_mean   = mean(Rsum_active,2);
Rsum_passive_mean  = mean(Rsum_passive,2);
Rsum_random_mean   = mean(Rsum_random,2);
Rsum_noRIS_mean    = mean(Rsum_noRIS,2);
Rsum_hybrid_iter_mean = mean(Rsum_hybrid_iter,2);

%% Plot
figure; hold on; box on; grid on;
plot(Ps_max_dBm,Rsum_active_mean,'-r^','LineWidth',1.5);
plot(Ps_max_dBm,Rsum_passive_mean,'-bo','LineWidth',1.5);
plot(Ps_max_dBm,Rsum_random_mean,'-m^','LineWidth',1.5);
plot(Ps_max_dBm,Rsum_noRIS_mean,'--k','LineWidth',1.5);
plot(Ps_max_dBm,Rsum_hybrid_iter_mean,'-g*','LineWidth',1.5);

xlabel('Total BS transmit power $P_{\rm s}^{\rm max}$ (dBm)','Interpreter','latex');
ylabel('Average sum-rate (bps/Hz)','Interpreter','latex');
title('');

labels = {'Active IRS (Optimized)','Passive IRS (Optimized)','Passive IRS (Random Phase)','Without IRS','Hybrid IRS (Optimized)'};
legend(labels,'Interpreter','latex','FontSize',12);

set(gca,'FontName','Times','FontSize',12);

save('sumrate_vs_power_with_hybrid_practicalTheta.mat','Ps_max_dBm','Ps_W',...
     'Rsum_active_mean','Rsum_passive_mean','Rsum_random_mean','Rsum_noRIS_mean', ...
     'Rsum_hybrid_iter_mean','Q');

toc
