clear

%% Set parameters

% Tensor size
d = 10;
n = 50 * ones(d, 1);

% number of tests
S = 7;

%% Tensor ranks setup and form the tensor

% Tensor rank
r_tmp_y = 5;
r_tmp_z = 5;

% generate temp tensor
tmp_y = TTrand(n, r_tmp_y);
tmp_z = TTrand(n, r_tmp_z);

% generate the tensors that used to do the hadamard product
% TT ranks of Y and Z are 25
% but thier minimal ranks are 5
% and A = Y \odot Z has a rank of 25*25=625
% but A's minimal rank is 5*5=25

y = tmp_y; z = tmp_z;
for i = 1 : 5 - 1
  y = y + tmp_y;
  z = z + tmp_z;
end

% ranks to test every 10
test_ranks = 5:5:60;

%% init matrix of errors and matrix of times

% we'll use the following functions to round A:
% 1. hadamard_round_randorth (HRO)
% 2. hadamard_round_randorth_without_svd (HRO_no_svd)
% 3. round (ori)
% 4. round_randorth (randorth)
% 5. round_orthrand (orthrand)
% 1 and 2 are algorithms we formulated, 3 is from Oseledets (2011) and 4 and 5 are from Daas (2023)

L = length(test_ranks);
errors_HRO = zeros(L, S);
errors_HRO_no_svd = zeros(L, S);
errors_ori = zeros(L, S);
errors_randorth = zeros(L, S);
errors_orthrand = zeros(L, S);

time_HRO = zeros(L, S);
time_HRO_no_svd = zeros(L, S);
time_ori = zeros(L, S);
time_randorth = zeros(L, S);
time_orthrand = zeros(L, S);

%% Real result of the tensor
% TT rank of tt is 25, we use Frobenius norm to calculate the error
% the norm function is from the TT-Toolbox: norm.m
tt = 5 * tmp_y .* (5 * tmp_z);

%% Test the algorithms

for i = 1 : L
  ell = test_ranks(i);
  for j = 1 : S
    HRO = tic;
    x = hadamard_round_randorth(y, z, ell);
    time_HRO(i, j) = toc(HRO);
    errors_HRO(i, j) = norm(tt - x) / norm(tt);

    HRO_no_svd = tic;
    x = hadamard_round_randorth_without_svd(y, z, ell);
    time_HRO_no_svd(i, j) = toc(HRO_no_svd);
    errors_HRO_no_svd(i, j) = norm(tt - x) / norm(tt);

    ori = tic;
    x = round(y .* z, 1e-5, ell);
    time_ori(i, j) = toc(ori);
    errors_ori(i, j) = norm(tt - x) / norm(tt);

    randorth = tic;
    x = round_randorth(y .* z, ell);
    time_randorth(i, j) = toc(randorth);
    errors_randorth(i, j) = norm(tt - x) / norm(tt);

    orthrand = tic;
    x = round_orthrand(y .* z, ell);
    time_orthrand(i, j) = toc(orthrand);
    errors_orthrand(i, j) = norm(tt - x) / norm(tt);
  end
end

%% Plot results

% Post-process errors

[err_HRO, neg_HRO, pos_HRO] = computeError(errors_HRO);
[err_HRO_no_svd, neg_HRO_no_svd, pos_HRO_no_svd] = computeError(errors_HRO_no_svd);
[err_ori, neg_ori, pos_ori] = computeError(errors_ori);
[err_randorth, neg_randorth, pos_randorth] = computeError(errors_randorth);
[err_orthrand, neg_orthrand, pos_orthrand] = computeError(errors_orthrand);

f = figure(1);
f.Position(1:2) = [0,1050];
f.Position(3:4) = [1050, 700];

subplot(1, 2, 1)

errorbar(test_ranks, err_ori, neg_ori, pos_ori,     'bo','markersize',10,'linewidth',2)
ax = gca;
ax.YScale = 'log';

hold on

errorbar(test_ranks, err_HRO, neg_HRO, pos_HRO,     'ks','markersize',10,'linewidth',2)
errorbar(test_ranks, err_HRO_no_svd, neg_HRO_no_svd, pos_HRO_no_svd, 'rx','markersize',10,'linewidth',2)
errorbar(test_ranks, err_randorth, neg_randorth, pos_randorth, 'g^','markersize',10,'linewidth',2)
errorbar(test_ranks, err_orthrand, neg_orthrand, pos_orthrand, 'm+','markersize',10,'linewidth',2)

hold off

% title('Relative error')
xlabel('Maximum Target Rank', 'FontSize', 18)
ylabel('Relative Error', 'FontSize', 18)
legend('ori', 'HRO', 'HRO-no-SVD', 'randorth', 'orthrand')
set(gca,'FontSize',16)

% Post-process timings

[speedup_HRO, neg_HRO, pos_HRO] = computeSpeedup(time_HRO, time_ori);
[speedup_HRO_no_svd, neg_HRO_no_svd, pos_HRO_no_svd] = computeSpeedup(time_HRO_no_svd, time_ori);
[speedup_ori, neg_ori, pos_ori] = computeSpeedup(time_ori, time_ori);
[speedup_randorth, neg_randorth, pos_randorth] = computeSpeedup(time_randorth, time_ori);
[speedup_orthrand, neg_orthrand, pos_orthrand] = computeSpeedup(time_orthrand, time_ori);

subplot(1, 2, 2)

errorbar(test_ranks, speedup_ori,       neg_ori,       pos_ori,       'bo-','markersize',10,'linewidth',2);

hold on

errorbar(test_ranks, speedup_HRO,       neg_HRO,       pos_HRO,       'ks-','markersize',10,'linewidth',2);
errorbar(test_ranks, speedup_HRO_no_svd, neg_HRO_no_svd, pos_HRO_no_svd, 'rx-','markersize',10,'linewidth',2);
errorbar(test_ranks, speedup_randorth,  neg_randorth,  pos_randorth,  'g^-','markersize',10,'linewidth',2);
errorbar(test_ranks, speedup_orthrand,  neg_orthrand,  pos_orthrand,  'm+-','markersize',10,'linewidth',2);

hold off
title("Observed speedup relative to TT-rounding")
xlabel('Maximum Target Rank', 'FontSize', 18)
ylabel('Speedup', 'FontSize', 18)
legend('ori', 'HRO', 'HRO-no-SVD', 'randorth', 'orthrand')
set(gca,'FontSize',16)

%% Custom function to post-process errors

function [error, neg, pos] = computeError(errors)
  % errors = errors(:,2:end);
  error = squeeze(median(errors, 2));
  neg = error -  squeeze(min(errors,[],2));
  pos = squeeze(max(errors,[],2)) - error;
end

%% Custom function to post-process runtimes

function [speedup, neg, pos] = computeSpeedup(time, timeref)
  % time = time(:,2:end,:);
  % timeref = timeref(:,2:end,:);
  speedup = median(timeref, 2) ./ median(time, 2);
  neg = speedup - median(timeref,2) ./ max(time,[],2);
  pos = median(timeref,2) ./ min(time,[],2) - speedup;
end