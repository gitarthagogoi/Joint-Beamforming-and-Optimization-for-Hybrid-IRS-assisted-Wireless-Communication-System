function [h_k, f_k, G, h_k_hat, f_k_hat, G_hat, sigma_e2] = Channel_generate2_imperfect(...
        K, N, M, large_fading_AI, large_fading_DI, Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser, f_c)
% CHANNEL_GENERATE2_IMPERFECT
% Generates true channels (h_k, f_k, G) and estimated versions (h_k_hat, f_k_hat, G_hat)
% using a controllable imperfect CSI model:
%
%    H_hat = accuracy * H + sqrt(1 - accuracy^2) * W,   W ~ CN(0,1)
%
% Outputs:
%  - h_k (K x M) true BS->User channels
%  - f_k (K x N) true RIS->User channels
%  - G   (N x M) true BS->RIS channel
%  - h_k_hat, f_k_hat, G_hat : estimated channels (imperfect CSI)
%  - sigma_e2 : struct with error variances (sigma_e2 = 1 - accuracy^2)

%% ----- Set estimation accuracy (edit as needed) -----
% Accuracy between 0 (completely wrong) and 1 (perfect)
accuracy_BSUser    = 0.80;   % BS -> User channel estimate accuracy
accuracy_RIStoUser = 0.70;   % RIS -> User channel estimate accuracy
accuracy_BStoRIS   = 0.90;   % BS -> RIS channel estimate accuracy
% -----------------------------------------------------

% Clamp to valid range
accuracy_BSUser    = min(max(accuracy_BSUser,0),1);
accuracy_RIStoUser = min(max(accuracy_RIStoUser,0),1);
accuracy_BStoRIS   = min(max(accuracy_BStoRIS,0),1);

% Derived estimation error variances
sigma_e2.BSUser    = 1 - accuracy_BSUser^2;
sigma_e2.RIStoUser = 1 - accuracy_RIStoUser^2;
sigma_e2.BStoRIS   = 1 - accuracy_BStoRIS^2;

%% ----- Allocate true-channel matrices -----
h_k = zeros(K,M);
f_k = zeros(K,N);
G   = zeros(N,M);

% === Generate true channels using your "channel_*2" definitions ===
for k = 1:K
    % BS -> User (use channel_H2 for this variant)
    h_k(k,:) = channel_H2(M,1,Dis_BStoUser(k),large_fading_AI,f_c);
end

for k = 1:K
    % RIS -> User
    f_k(k,:) = channel_F(N,1,Dis_RIStoUser(k),large_fading_DI,f_c);
end

% BS -> RIS
G(:,:) = channel_G(N,M,Dis_BStoRIS,large_fading_DI,f_c);

%% ----- Create estimated channels (imperfect CSI) -----
% Complex Gaussian noise matrices for estimation error
W_h = (randn(K,M) + 1j*randn(K,M))/sqrt(2);
W_f = (randn(K,N) + 1j*randn(K,N))/sqrt(2);
W_G = (randn(N,M) + 1j*randn(N,M))/sqrt(2);

% Estimated channels
h_k_hat = accuracy_BSUser    .* h_k + sqrt(sigma_e2.BSUser)    .* W_h;
f_k_hat = accuracy_RIStoUser .* f_k + sqrt(sigma_e2.RIStoUser) .* W_f;
G_hat   = accuracy_BStoRIS   .* G   + sqrt(sigma_e2.BStoRIS)   .* W_G;

end
