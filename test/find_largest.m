clear
ell = 5;
D_set = 10:10:50;
S = 7;
L = length(D_set);
for fun_num = 1:2
  if fun_num == 1
    t = linspace(-500, 500, 10);
  else
    t = linspace(-2.5*pi, 2.5*pi, 10);
  end
  errors_HRO = zeros(L, S);
  errors_HRO_no_svd = zeros(L, S);
  errors_TTrounding = zeros(L, S);
  errors_randorth = zeros(L, S);
  errors_orthrand = zeros(L, S);
  errors_twosided = zeros(L, S);
  
  time_HRO = zeros(L, S);
  time_HPCRL = zeros(L, S);
  time_HBilinearForm = zeros(L, S);
  
  time_HRO_no_svd = zeros(L, S);
  time_HPCRL_no_svd = zeros(L, S);
  time_HBilinearForm_no_svd = zeros(L, S);
  
  time_TTrounding = zeros(L, S);
  
  time_randorth = zeros(L, S);
  time_PCRL = zeros(L, S);
  time_dot = zeros(L, S);
  
  time_orthrand = zeros(L, S);
  
  time_twosided = zeros(L, S);
  
  for i = 1 : L
    D = D_set(i);
    if fun_num == 1
      A = qing(t, D);
    else
      A = alpine(t, D);
    end
    x = TTrand(length(t) * ones(1, D), 1);
    x.core = ones(1, x.ps(end) - 1);
    if fun_num == 1
      act_M = sum((500 ^ 2 - (1:D)) .^ 2);
    else
      act_M = D * 2.6*pi;
    end
    x1 = x;
    iter_num = 100;
    %% start iteration
    for j = 1:S
      HRO = tic;
      for iter = 1:iter_num
        [x, time_HPCRL(i, j)] = hadamard_round_randorth(A, x, ell);
        x = x/norm(x);
      end
      HBF = tic;
      M_svd = HBilinearForm(A, x);
      time_HBilinearForm(i, j) = toc(HBF);
      time_HRO(i, j) = toc(HRO);
      errors_HRO(i, j) = abs(M_svd - act_M)/act_M ;
  
      x = x1;
      HRO_no_svd = tic;
      for iter = 1:iter_num
        [x, time_HPCRL_no_svd(i, j)] = hadamard_round_randorth_without_svd(A, x, ell);
        x = x/norm(x);
      end
      HBF_no_svd = tic;
      M_no_svd = HBilinearForm_without_svd(A, x);
      time_HBilinearForm_no_svd(i, j) = toc(HBF_no_svd);
      time_HRO_no_svd(i, j) = toc(HRO_no_svd);
      errors_HRO_no_svd(i, j) = abs(M_no_svd - act_M)/act_M;
  
      x = x1;
      randorth = tic;
      for iter = 1:iter_num
        [x, time_PCRL(i, j)] = round_randorth(A.*x, ell);
        x = x/norm(x);
      end
      DOT = tic;
      M_randorth = dot(x, A.*x);
      time_dot(i, j) = toc(DOT);
      time_randorth(i, j) = toc(randorth);
      errors_randorth(i, j) = abs(M_randorth - act_M)/act_M;
  
      x = x1;
      TTRounding = tic;
      for iter = 1:iter_num
        x = round(A.*x, eps, ell);
        x = x/norm(x);
      end
      M_TTRounding = dot(x, A.*x);
      time_TTrounding(i, j) = toc(TTRounding);
      errors_TTrounding(i, j) = abs(M_TTRounding - act_M)/act_M;
  
      x = x1;
      orthrand = tic;
      for iter = 1:iter_num
        x = round_orthrand(A.*x, ell);
        x = x/norm(x);
      end
      M_orthrand = dot(x, A.*x);
      time_orthrand(i, j) = toc(orthrand);
      errors_orthrand(i, j) = abs(M_orthrand - act_M)/act_M;
  
      x = x1;
      twosided = tic;
      for iter = 1:iter_num
        x = round_twosided(A.*x, ell);
        x = x/norm(x);
      end
      M_twosided = dot(x, A.*x);
      time_twosided(i, j) = toc(twosided);
      errors_twosided(i, j) = abs(M_twosided - act_M)/act_M;
    end
  end
  if fun_num == 1
    save('qing.mat')
  else
    save('alpine.mat')
  end
end

function y = qing(t, D)
  tmp = tt_tensor();
  n = length(t) * ones(1, D);
  r = [1, ones(1, D - 1), 1];
  ps = cumsum([1, n .* r(1 : end - 1) .* r(2 : end)]);
  cr = zeros(ps(end) - 1, 1);
  tmp.n = n;
  tmp.d = D;
  tmp.r = r;
  tmp.ps = ps;
  y.core = cr;
  cr = cr + 1;
  tmp.core = cr;
  for k = 1 : D
    tmp.core(ps(k) : ps(k + 1) - 1) = (t .^ 2 - k) .^ 2;
    if k == 1
        y = tmp;
    else
        y = y + tmp;
    end
    tmp.core = cr;
  end
end

function y = alpine(t, D)
  tmp = tt_tensor();
  n = length(t) * ones(1, D);
  r = [1, ones(1, D - 1), 1];
  ps = cumsum([1, n .* r(1 : end - 1) .* r(2 : end)]);
  cr = zeros(ps(end) - 1, 1);
  tmp.n = n;
  tmp.d = D;
  tmp.r = r;
  tmp.ps = ps;
  y.core = cr;
  cr = cr + 1;
  tmp.core = cr;
  for k = 1 : D
    tmp.core(ps(k) : ps(k + 1) - 1) = abs(t .* sin(t) + 0.1 * t);
    if k == 1
        y = tmp;
    else
        y = y + tmp;
    end
    tmp.core = cr;
  end
end