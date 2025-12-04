function [Rsum_DF, stats] = DF_precoding_full(K, Ps_max, Pr_max, sigma2, h_k, f_k, G, varargin)
% Decode-and-Forward (DF) relaying sum-rate following Björnson et al. Eq. (10)
%
% Inputs:
%   K         - Number of users
%   Ps_max    - Source (BS) average transmit power (Watts)
%   Pr_max    - Relay (RIS) power budget (Watts)
%   sigma2    - Noise power (Watts)
%   h_k       - BS->User channel matrix (K x M)
%   f_k       - RIS->User channel matrix (K x N)
%   G         - BS->RIS channel matrix (N x M)
%
% Optional flags (Name,Value):
%   'optPower' (true/false) : Enable optimal p1/p2 from Proposition 1
%   'clipSNR'  (true/false) : Clip extreme SNRs for realism
%
% Outputs:
%   Rsum_DF - total DF sum-rate (bps/Hz)
%   stats   - structure with per-user SNRs and power allocation

%% Parse optional arguments
p = inputParser;
addParameter(p, 'optPower', true, @islogical);
addParameter(p, 'clipSNR', true, @islogical);
parse(p, varargin{:});
optPower = p.Results.optPower;
clipSNR  = p.Results.clipSNR;

%% Dimensions & preparation
[N, M] = size(G);
[Kh, Mh] = size(h_k);
[Kf, Nf] = size(f_k);
if Mh == 1 && M > 1, h_k = repmat(h_k, 1, M); end
if Nf == 1 && N > 1, f_k = repmat(f_k, 1, N); end

Rsum_DF = 0;
stats.gamma_SR = zeros(K,1);
stats.gamma_SD = zeros(K,1);
stats.gamma_RD = zeros(K,1);
stats.p1 = zeros(K,1);
stats.p2 = zeros(K,1);

%% Main per-user DF computation
for k = 1:K
    % --- Channel average gains (? parameters) ---
    beta_sr = norm(G,'fro')^2 / (M*N);            % BS->Relay
    beta_rd = norm(f_k(k,:))^2 / N;               % Relay->User
    beta_sd = norm(h_k(k,:))^2 / max(1,size(h_k,2)); % Direct BS->User

    % --- Power allocation (Proposition 1) ---
    if optPower
        p_avg = Ps_max;  % interpret as average power
        if beta_sd > beta_sr
            % Relay suboptimal ? use SISO mode
            p1 = 2*p_avg; p2 = 0;
        else
            denom = (beta_sr + beta_rd - beta_sd);
            if denom <= 0
                p1 = 2*p_avg; p2 = 0;
            else
                p1 = 2*p_avg * beta_rd / denom;
                p2 = 2*p_avg * (beta_sr - beta_sd) / denom;
            end
            p1 = max(p1,0); p2 = max(p2,0);
        end
    else
        p1 = Ps_max;
        p2 = Pr_max;
    end

    % --- Compute SNRs ---
    gamma_SR = p1 * beta_sr / sigma2;      % first hop
    gamma_SD = p1 * beta_sd / sigma2;      % direct link
    gamma_RD = p2 * beta_rd / sigma2;      % relay link
    gamma_2hop = gamma_SD + gamma_RD;      % combined second hop

    % --- Optional clipping (avoid unrealistically huge SNRs) ---
    if clipSNR
        maxSNR = 1e4;   % cap at ~40 dB
        gamma_SR = min(gamma_SR, maxSNR);
        gamma_2hop = min(gamma_2hop, maxSNR);
    end

    % --- DF rate (half-duplex factor 0.5) ---
    Rk = 0.5 * log2(1 + min(gamma_SR, gamma_2hop));
    Rsum_DF = Rsum_DF + Rk;

    % Store debug stats
    stats.gamma_SR(k) = gamma_SR;
    stats.gamma_SD(k) = gamma_SD;
    stats.gamma_RD(k) = gamma_RD;
    stats.p1(k) = p1;
    stats.p2(k) = p2;
end
end
