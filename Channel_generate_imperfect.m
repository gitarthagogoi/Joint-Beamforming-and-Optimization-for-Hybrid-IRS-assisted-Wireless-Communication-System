function [h_k, f_k, G, h_k_hat, f_k_hat, G_hat, sigma_e2] = Channel_generate_imperfect(...
        K, N, M, large_fading_AI, large_fading_DI, Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser, f_c)
% CHANNEL_GENERATE_IMPERFECT (no opts input)
% Generates true channels and estimated channels using a static, in-function
% "accuracy" mixing model:
%    H_hat = accuracy * H + sqrt(1 - accuracy^2) * W,   W ~ CN(0,1)
%
% Edit the accuracy variables below to change estimation quality.
%
% Outputs:
%  - h_k (K x M) true BS->User channels
%  - f_k (K x N) true RIS->User channels
%  - G   (N x M) true BS->RIS channel
%  - h_k_hat, f_k_hat, G_hat : estimated channels
%  - sigma_e2 : struct with error variances (sigma_e2 = 1 - accuracy^2)

%% ----- Set estimation accuracy here (edit these values) -----
% You can either set a single global accuracy or per-link accuracies.
% Valid range: 0 (no information) to 1 (perfect CSI).

% Option A: single global accuracy (uncomment to use)
% global_accuracy = 0.80;   % 80% accurate for all links

% Option B: per-link accuracies (preferred if you want different values)
accuracy_BSUser    = 0.80;   % BS -> User channel estimate accuracy (0..1)
accuracy_RIStoUser = 0.70;   % RIS -> User channel estimate accuracy (0..1)
accuracy_BStoRIS   = 0.90;   % BS -> RIS channel estimate accuracy (0..1)

% If you want to use global_accuracy instead of per-link, uncomment these lines:
% accuracy_BSUser = global_accuracy;
% accuracy_RIStoUser = global_accuracy;
% accuracy_BStoRIS = global_accuracy;
%% ---------------------------------------------------------------

% clamp accuracies to [0,1] just in case
accuracy_BSUser    = min(max(accuracy_BSUser,0),1);
accuracy_RIStoUser = min(max(accuracy_RIStoUser,0),1);
accuracy_BStoRIS   = min(max(accuracy_BStoRIS,0),1);

% derived error variances (for convenience)
sigma_e2.BSUser    = 1 - accuracy_BSUser^2;
sigma_e2.RIStoUser = 1 - accuracy_RIStoUser^2;
sigma_e2.BStoRIS   = 1 - accuracy_BStoRIS^2;

%% ----- allocate true-channel matrices -----
h_k = zeros(K,M);   % each row: user k, columns: M BS antennas
f_k = zeros(K,N);   % each row: user k, columns: N RIS elements
G   = zeros(N,M);   % N x M BS->RIS

% generate true channels using your helper functions
for k = 1:K
    h_k(k,:) = channel_H(M,1,Dis_BStoUser(k),large_fading_AI,f_c);
end

for k = 1:K
    f_k(k,:) = channel_F(N,1,Dis_RIStoUser(k),large_fading_DI,f_c);
end

G(:,:) = channel_G(N,M,Dis_BStoRIS,large_fading_DI,f_c);

%% ----- create estimated channels using the accuracy mixing model -----
% Complex Gaussian noise W of same size as each channel
W_h = (randn(K,M) + 1j*randn(K,M))/sqrt(2);
W_f = (randn(K,N) + 1j*randn(K,N))/sqrt(2);
W_G = (randn(N,M) + 1j*randn(N,M))/sqrt(2);

% estimated channels
h_k_hat = accuracy_BSUser    .* h_k + sqrt(sigma_e2.BSUser)    .* W_h;
f_k_hat = accuracy_RIStoUser .* f_k + sqrt(sigma_e2.RIStoUser) .* W_f;
G_hat   = accuracy_BStoRIS   .* G   + sqrt(sigma_e2.BStoRIS)   .* W_G;

end
