clear

%% Set parameters

% Tensor size
d = 10;
n = 50 * ones(1, d);

% number of tests
S = 6;

%% Tensor ranks setup and form the tensor

% Tensor rank
r_tmp_y = 5;
r_tmp_z = 5;

% generate temp tensor
tmp_y = TTrandn(n, r_tmp_y);
tmp_z = TTrandn(n, r_tmp_z);

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

% ranks to test every 5
test_ranks = 30:10:120;

%% init matrix of errors and matrix of times

% we'll use the following functions to round A:
% 1. hadamard_round_randorth (HaTT)
% 2. hadamard_round_randorth_without_svd (HaTT2)
% 3. round (TTrounding)
% 4. round_randorth (randorth)
% 5. round_orthrand (orthrand)
% 1 and 2 are algorithms we formulated, 3 is from Oseledets (2011) and 4 and 5 are from Daas (2023)

L = length(test_ranks);
errors_HaTT1 = zeros(L, S);
errors_HaTT2 = zeros(L, S);
% errors_TTrounding = zeros(L, S);
% errors_randorth = zeros(L, S);
% errors_orthrand = zeros(L, S);

time_HaTT1 = zeros(L, S);
time_HaTT2 = zeros(L, S);
time_HPCRL = zeros(L, S);
time_HPCRL_no_svd = zeros(L, S);
% time_TTrounding = zeros(L, S);
% time_randorth = zeros(L, S);
% time_orthrand = zeros(L, S);

%% Real result of the tensor
% TT rank of tt is 25, we use Frobenius norm to calculate the error
% the norm function is from the TT-Toolbox: norm.m
tt = (5*tmp_y) .* (5*tmp_z);

%% Test the algorithms

for i = 1 : L
  ell = test_ranks(i);
  for j = 1 : S
    HaTT1 = tic;
    [x, time_HPCRL(i, j)] = hadamard_round_randorth(y, z, ell);
    time_HaTT1(i, j) = toc(HaTT);
    errors_HaTT1(i, j) = norm(tt - x) / norm(tt);

    HaTT2 = tic;
    [x, time_HPCRL_no_svd(i, j)] = hadamard_round_randorth_without_svd(y, z, ell);
    time_HaTT2(i, j) = toc(HaTT2);
    errors_HaTT2(i, j) = norm(tt - x) / norm(tt);

    % TTrounding = tic;
    % x = round(y .* z, 1e-5, ell);
    % time_TTrounding(i, j) = toc(TTrounding);
    % errors_TTrounding(i, j) = norm(tt - x) / norm(tt);
    % 
    % randorth = tic;
    % x = round_randorth(y .* z, ell);
    % time_randorth(i, j) = toc(randorth);
    % errors_randorth(i, j) = norm(tt - x) / norm(tt);
    % 
    % orthrand = tic;
    % x = round_orthrand(y .* z, ell);
    % time_orthrand(i, j) = toc(orthrand);
    % errors_orthrand(i, j) = norm(tt - x) / norm(tt);
  end
end

%% save data for reusing and plot result
save("randtest.mat");

% PlotResult("randtest.mat")