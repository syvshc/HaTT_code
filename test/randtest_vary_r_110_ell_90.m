% clear

%% Set parameters
load('randtest_vary_r_ell_90.mat')
% Tensor size
d = 7;
n = 20 * ones(d, 1);

% number of tests
S = 7;

%% Tensor ranks setup and form the tensor

% test tensor ranks every 5
test_ranks = 60:10:150;

% generate temp tensor

% generate the tensors that used to do the hadamard product
% TT ranks of Y and Z are 25
% but thier minimal ranks are 5
% and A = Y \odot Z has a rank of 25*25=625
% but A's minimal rank is 5*5=25


% we fix the target rank 50
ell = 90;

%% init matrix of errors and matrix of times

% we'll use the following functions to round A:
% 1. round (TTrounding)
% 2. round_randorth (randorth)
% 3. round_orthrand (orthrand)
% 4. round_twosided (twosided)
% 5. HaTT1 (HaTT)
% 6. HaTT2 (HaTT_no_svd)
% 5 and 6 are algorithms we formulated, 1 is from Oseledets (2011) and 2, 3 and 4 are from Daas (2023)

L = length(test_ranks);
% errors_HaTT1 = zeros(L, S);
% errors_HaTT2 = zeros(L, S);
% errors_TTrounding = zeros(L, S);
% errors_randorth = zeros(L, S);
% errors_orthrand = zeros(L, S);
% errors_twosided = zeros(L, S);
% time_HaTT1 = zeros(L, S);
% time_HaTT2 = zeros(L, S);
% time_TTrounding = zeros(L, S);
% time_randorth = zeros(L, S);
% time_orthrand = zeros(L, S);
% time_HPCRL = zeros(L, S);
% time_HPCRL_no_svd = zeros(L, S);
% time_PCRL = zeros(L, S);
% time_HMcore0 = zeros(L, S);
% time_HMcore0_no_svd = zeros(L, S);
% time_Mcore0 = zeros(L, S);
% time_hadamard = zeros(L, 1);
% time_twosided = zeros(L, S);

%% Test the algorithms
i=6;
% for i = 1 : L
  r = test_ranks(i);
  y = TTrand(n, r);
  z = TTrand(n, r);
%   try
%     hadamard = tic;
%     tt = y .* z;
%     time_hadamard(i) = toc(hadamard);
%     normtt = norm(tt);
%   catch
%     r_overflow = r;
%     try
%       tt = [0; 1];
%       TT = full(y) .* full(z);
%       normTT = norm(TT);
%     catch
%       TT = [0; 1];
%     end
%   end
  % only for r = 110
  try
%       tt = y .* z;
      TT = full(y) .* full(z);
      normTT = norm(TT);
  catch
      TT = [0; 1];
  end
  for j = 1 : S
    try
      TTrounding = tic;
      x = round(tt, 1e-10, ell);
      time_TTrounding(i, j) = toc(TTrounding);
    catch
      time_TTrounding(i, j) = inf;
      x = [0; 1; 1];
    end
%     try
%       errors_TTrounding(i, j) = norm(tt - x) / normtt;
%     catch
      try
        X = full(x);
        errors_TTrounding(i, j) = norm(TT-X) / normTT;
      catch
        errors_TTrounding(i, j) = inf;
      end
%     end

    try
      HaTT = tic;
      [x, time_HPCRL(i, j), time_HMcore0(i, j)] = HaTT1(y, z, ell);
      time_HaTT1(i, j) = toc(HaTT);
    catch
      time_HaTT1(i, j) = inf;
      x = [0; 1; 1];
    end
%     try
%       errors_HaTT1(i, j) = norm(tt - x) / normtt;
%     catch
      try
        X = full(x);
        errors_HaTT1(i, j) = norm(TT-X) / normTT;
      catch
        errors_HaTT1(i, j) = inf;
      end
%     end

    try
      HaTT_no_svd = tic;
      [x, time_HPCRL_no_svd(i, j), time_HMcore0_no_svd(i, j)] = HaTT2(y, z, ell);
      time_HaTT2(i, j) = toc(HaTT_no_svd);
    catch
      time_HaTT2(i, j) = inf;
      x = [0; 1; 1];
    end
%     try
%       errors_HaTT2(i, j) = norm(tt - x) / normtt;
%     catch
      try
        X = full(x);
        errors_HaTT2(i, j) = norm(TT-X) / normTT;
      catch
        errors_HaTT2(i, j) = inf;
      end
%     end
    
    try
      randorth = tic;
      [x, time_PCRL(i, j), time_Mcore0(i, j)] = round_randorth(tt, ell);
      time_randorth(i, j) = toc(randorth);
    catch
      time_randorth(i, j) = inf;
      x = [0; 1; 1];
    end
%     try
%       errors_randorth(i, j) = norm(tt - x) / normtt;
%     catch
      try
        X = full(x);
        errors_randorth(i, j) = norm(TT-X) / normTT;
      catch
        errors_randorth(i, j) = inf;
      end
%     end
    
    try
      orthrand = tic;
      x = round_orthrand(tt, ell);
      time_orthrand(i, j) = toc(orthrand);
    catch
      time_orthrand(i, j) = inf;
      x = [0; 1; 1];
    end
%     try
%       errors_orthrand(i, j) = norm(tt - x) / normtt;
%     catch
      try
        X = full(x);
        errors_orthrand(i, j) = norm(TT-X) / normTT;
      catch
        errors_orthrand(i, j) = inf;
      end
%     end

    try
      twosided = tic;
      x = round_twosided(tt, ell);
      time_twosided(i, j) = toc(twosided);
    catch
      time_twosided(i, j) = inf;
      x = [0; 1; 1];
    end
%     try
%       errors_twosided(i, j) = norm(tt - x) / normtt;
%     catch
      try
        X = full(x);
        errors_twosided(i, j) = norm(TT-X) / normTT;
      catch
        errors_twosided(i, j) = inf;
      end
%     end
  end
% end

%% save data for reusing and plot result
save("randtest_vary_r_ell_90.mat");
% 
PlotResult_vary_r("randtest_vary_r_ell_90.mat")