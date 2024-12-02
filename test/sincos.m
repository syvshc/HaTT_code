clear

item_num = 60;
d = 7;
n=15;
len = n^d;
t1 = 2*pi/len:2*pi/len:2*pi;
a = 10 * rand(1, item_num) + 0.1;
b = 10 * rand(1, item_num) + 0.1;
y1 = 0; z1 = 0;
for k = 1: item_num
  y1 = y1 + a(k) * sin(k * t1);
  z1 = z1 + b(k) * cos(k * t1);
end
y = y1; z = z1;
clear y1 z1 t1 a b;

sz = n * ones(1, d);
y = tt_tensor(reshape(y, sz))
z = tt_tensor(reshape(z, sz))

% ranks to test every 2
test_ranks = 4:4:60;
S = 7;
%% init matrix of errors and matrix of times

% we'll use the following functions to round A:
% 1. round (TTrounding)
% 2. round_randorth (randorth)
% 3. round_orthrand (orthrand)
% 4. round_twosided (twosided)
% 5. HaTT1 (HaTT)
% 6. HaTT2 (HaTT_no_svd)
% 5 and 6 are algorithms we formulated, 1 is from Oseledets (2011) and 2, 3, 4 are from Daas (2023)

L = length(test_ranks);
errors_HaTT1 = zeros(L, S);
errors_HaTT2 = zeros(L, S);
errors_TTrounding = zeros(L, S);
errors_randorth = zeros(L, S);
errors_orthrand = zeros(L, S);
errors_twosided = zeros(L, S);
time_HaTT1 = zeros(L, S);
time_HaTT2 = zeros(L, S);
time_TTrounding = zeros(L, S);
time_randorth = zeros(L, S);
time_orthrand = zeros(L, S);
time_twosided = zeros(L, S);

%% Real result of the tensor
% the norm function is from the TT-Toolbox: norm.m
tt = y .* z

%% Test the algorithms

for i = 1 : L
  ell = test_ranks(i);
  for j = 1 : S
    HaTT = tic;
    x = HaTT1(y, z, ell);
    time_HaTT1(i, j) = toc(HaTT);
    errors_HaTT1(i, j) = norm(tt - x) / norm(tt);

    HaTT_no_svd = tic;
    x = HaTT2(y, z, ell);
    time_HaTT2(i, j) = toc(HaTT_no_svd);
    errors_HaTT2(i, j) = norm(tt - x) / norm(tt);

    TTrounding = tic;
    x = round(y .* z, 1e-16, ell);
    time_TTrounding(i, j) = toc(TTrounding);
    errors_TTrounding(i, j) = norm(tt - x) / norm(tt);

    randorth = tic;
    x = round_randorth(y .* z, ell);
    time_randorth(i, j) = toc(randorth);
    errors_randorth(i, j) = norm(tt - x) / norm(tt);

    orthrand = tic;
    x = round_orthrand(y .* z, ell);
    time_orthrand(i, j) = toc(orthrand);
    errors_orthrand(i, j) = norm(tt - x) / norm(tt);

    twosided = tic;
    x = round_twosided(y .* z, ell);
    time_twosided(i, j) = toc(twosided);
    errors_twosided(i, j) = norm(tt - x) / norm(tt);
  end
end

%% save data for reusing and plot result
save("tmp.mat");

PlotResult_sincos("tmp.mat")
