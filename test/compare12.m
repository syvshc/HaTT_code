clear

%% Set parameters

% We will format a Hilbert-type TT tensor

n = 20 * ones(7, 1);
d = length(n);
r = [1; 100 * ones(d - 1, 1); 1];

%% Set cores

y = tt_tensor();
y.n = n; y.d = d; y.r = r;
y.ps = cumsum([1; r(1 : end-1) .* n .* r(2 : end)]);
y.core = zeros(y.ps(end)-1, 1);
for k = 1 : d
  core = zeros(r(k), n(k), r(k + 1));
  for i1 = 1 : r(k)
    for i2 = 1 : n(k)
      for i3 = 1 : r(k + 1)
        core(i1, i2, i3) = 1 / (i1 + i2 + i3 - 1);
      end
    end
  end
  y.core(y.ps(k) : y.ps(k + 1) - 1) = core(:);
end
z = y;
% Tensor size
% d = 10;
% n = 50 * ones(1, d);

% number of tests
S = 7;

% %% Tensor ranks setup and form the tensor

% % Tensor rank
% r_tmp_y = 5;
% r_tmp_z = 5;

% % generate temp tensor
% tmp_y = TTrandn(n, r_tmp_y);
% tmp_z = TTrandn(n, r_tmp_z);

% % generate the tensors that used to do the hadamard product
% % TT ranks of Y and Z are 25
% % but thier minimal ranks are 5
% % and A = Y \odot Z has a rank of 25*25=625
% % but A's minimal rank is 5*5=25

% y = tmp_y; z = tmp_z;
% for i = 1 : 5 - 1
%   y = y + tmp_y;
%   z = z + tmp_z;
% end

% ranks to test every 5
trunc = 5;
test_ranks = 30:10:120;

%% init matrix of errors and matrix of times

% we'll use the following functions to round A:
% 1. HaTT1 (HaTT)
% 2. HaTT2 (HaTT_no_svd)
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
time_HPCRL1 = zeros(L, S);
time_HPCRL2 = zeros(L, S);
% time_TTrounding = zeros(L, S);
% time_randorth = zeros(L, S);
% time_orthrand = zeros(L, S);

%% Real result of the tensor
% TT rank of tt is 25, we use Frobenius norm to calculate the error
% the norm function is from the TT-Toolbox: norm.m
% tt = y .* z;
% normtt = norm(tt);
Y = full(y);
TT = Y .* Y;
normTT = norm(TT);

%% Test the algorithms

for i = 1 : L
    i
  ell = test_ranks(i);
  for j = 1 : S
    HaTT = tic;
    [x, time_HPCRL1(i, j)] = HaTT1(y, z, ell, 'trunc', trunc);
    time_HaTT1(i, j) = toc(HaTT);
    X = full(x);
    errors_HaTT1(i, j) = norm(TT - X) / normTT;
%     errors_HaTT1(i, j) = norm(tt - x) / normtt;

    HaTT_no_svd = tic;
    [x, time_HPCRL2(i, j)] = HaTT2(y, z, ell);
    time_HaTT2(i, j) = toc(HaTT_no_svd);
    X = full(x);
    errors_HaTT2(i, j) = norm(TT - X) / normTT;
%     errors_HaTT1(i, j) = norm(tt - x) / normtt;
  end
end

%% save data for reusing and plot result
save("compare12.mat");

PlotResult_compare12("compare12.mat")