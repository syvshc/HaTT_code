clear
d_set = 5:6;

%% setup and initialize
% number of tests
S = 2;
L = length(d_set);

errors_randorth = zeros(L, S);
errors_orthrand = zeros(L, S);
errors_twosided = zeros(L, S);
errors_HaTT1 = zeros(L, S);
errors_HaTT2 = zeros(L, S);
time_TTrounding = zeros(L, S);
time_randorth = zeros(L, S);
time_orthrand = zeros(L, S);
time_twosided = zeros(L, S);
time_HaTT1 = zeros(L, S);
time_HaTT2 = zeros(L, S);
round_time_TTrounding = zeros(L, S);
round_time_randorth = zeros(L, S);
round_time_orthrand = zeros(L, S);
round_time_twosided = zeros(L, S);
round_time_HaTT1 = zeros(L, S);
round_time_HaTT2 = zeros(L, S);

for i = 1 : L
  d = d_set(i);
  n = 2^d;
  hx = 2 * pi / n;
  xx = hx * ((1 : n)' - 0.5);
  try 
    [X, Y, Z] = meshgrid(xx, xx, xx);
    ori_phi = 0.2 * sin(X) .* sin(Y) .* sin(Z);
    disp(['====d = ', num2str(d), ' is being tested.====']);
  catch
    disp(['!!!d = ', num2str(d), ' is too large.!!!']);
    break;
  end
  try 
    [round_time, time, phi] = qttRDVDAC3(d, ori_phi, "TTrounding");
    time_TTrounding(i, :) = time;
    round_time_TTrounding(i, :) = sum(round_time);
    phi_TTrounding = phi;
  catch
    time_TTrounding(i, :) = inf;
    round_time_TTrounding(i, :) = inf;
  end
  for j = 1:S
    disp(['S = ', num2str(j), ' is being tested.']);
    try
      [round_time, time, phi] = qttRDVDAC3(d, ori_phi, "randorth");
      time_randorth(i, j) = time;
      round_time_randorth(i, j) = sum(round_time);
      try
        errors_randorth(i, j) = norm(phi - phi_TTrounding) / norm(phi_TTrounding);
      catch
        errors_randorth(i, j) = -1;
      end
    catch
      time_randorth(i, j) = inf;
      round_time_randorth(i, j) = inf;
    end
    try
      [round_time, time, phi] = qttRDVDAC3(d, ori_phi, "orthrand");
      time_orthrand(i, j) = time;
      round_time_orthrand(i, j) = sum(round_time);
      try
        errors_orthrand(i, j) = norm(phi - phi_TTrounding) / norm(phi_TTrounding);
      catch
        errors_orthrand(i, j) = -1;
      end
    catch
      time_orthrand(i, j) = inf;
      round_time_orthrand(i, j) = inf;
    end
    try
      [round_time, time, phi] = qttRDVDAC3(d, ori_phi, "twosided");
      time_twosided(i, j) = time;
      round_time_twosided(i, j) = sum(round_time);
      try
        errors_twosided(i, j) = norm(phi - phi_TTrounding) / norm(phi_TTrounding);
      catch
        errors_twosided(i, j) = -1;
      end
    catch
      time_twosided(i, j) = inf;
      round_time_twosided(i, j) = inf;
    end
    try
      [round_time, time, phi] = HaTT1_qttRDVDAC3(d, ori_phi);
      time_HaTT1(i, j) = time;
      round_time_HaTT1(i, j) = sum(round_time);
      try
        errors_HaTT1(i, j) = norm(phi - phi_TTrounding) / norm(phi_TTrounding);
      catch
        errors_HaTT1(i, j) = -1;
      end
    catch
      time_HaTT1(i, j) = inf;
      round_time_HaTT1(i, j) = inf;
    end
    try
      [round_time, time, phi] = HaTT2_qttRDVDAC3(d, ori_phi);
      time_HaTT2(i, j) = time;
      round_time_HaTT2(i, j) = sum(round_time);
      try
        errors_HaTT2(i, j) = norm(phi - phi_TTrounding) / norm(phi_TTrounding);
      catch
        errors_HaTT2(i, j) = -1;
      end
    catch
      time_HaTT2(i, j) = inf;
      round_time_HaTT2(i, j) = inf;
    end
  end
end
save(qttRDVDAC3.mat, "-regexp", "^errors_", "^time_", "^round_time_");