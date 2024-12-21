clear
d_set = [8,9:3:21];

%% setup and initialize
% number of tests
S = 7;
L = length(d_set);
errors_TTrounding = zeros(L, S);
errors_randorth = zeros(L, S);
errors_orthrand = zeros(L, S);
errors_twosided = zeros(L, S);
errors_HaTT1 = zeros(L, S);
errors_HaTT2 = zeros(L, S);
time_pcg = zeros(L, S);
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
energy_pcg = cell(L, 1);
energy_TTrounding = cell(L, 1);
energy_randorth = cell(L, S);
energy_orthrand = cell(L, S);
energy_twosided = cell(L, S);
energy_HaTT1 = cell(L, S);
energy_HaTT2 = cell(L, S);

bool_energy = 1;
bool_save = 0;

for i = 1 : L
  d = d_set(i);
  n = 2^d;
  hx = 2 * pi / n;
  xx = hx * ((1 : n)' - 0.5);
  % generate the initial value (matrix form)
  try 
    ori_phi_mat = zeros(n, n, n);
    for X = 1:n
        for Y = 1:n
            for Z = 1:n
                ori_phi_mat(X, Y, Z) = 0.2 * sin(xx(X)) .* sin(xx(Y)) .* sin(xx(Z));
            end
        end
    end
    disp(['====d = ', num2str(d), ' is being tested.====']);
      [time, phi_pcg, energy] = RDVDAC3(d, ori_phi_mat, bool_energy, bool_save);
      time_pcg(i, :) = time;
      energy_pcg{i} = energy;
      phi_pcg = reshape(phi_pcg, 2 * ones(1, 3 * d));
      phi_pcg = tt_tensor(phi_pcg);
  catch
    disp(['!!!d = ', num2str(d), ' is too large for direct representation.!!!']);
    phi_pcg = [0, 1];
  end
  clearvars xx;
  % generate the initial value (TT form)
cores = cell(d, 1);
cores{1}(:, 1) = sin(hx*[0.5;1.5]);
cores{1}(:, 2) = cos(hx*[0.5;1.5]);
for j = 2 : (d - 1)
    cores{j}(:, 1, 1) = cos(2^(j-1)*hx*[0;1]);
    cores{j}(:, 2, 2) = cores{j}(:, 1, 1);
    cores{j}(:, 1, 2) = -sin(2^(j-1)*hx*[0;1]);
    cores{j}(:, 2, 1) = -cores{j}(:, 1, 2);
end
cores{end}(:, 1) = cos(2^(d-1)*hx*[0;1]);
cores{end}(:, 2) = sin(2^(d-1)*hx*[0;1]);
tt = tt_tensor(cores);
ori_phi_tt = tt_tensor;
ori_phi_tt.core = [tt.core;tt.core;tt.core];
ori_phi_tt.r = [tt.r;tt.r(2:end);tt.r(2:end)];
ori_phi_tt.n = [tt.n;tt.n;tt.n];
ori_phi_tt.d = length(ori_phi_tt.n);
ori_phi_tt.ps = cumsum([1;ori_phi_tt.n.*ori_phi_tt.r(1:end - 1).*ori_phi_tt.r(2:end)]);
ori_phi_tt = 0.2 * ori_phi_tt;
  try 
    [round_time, time, phi, energy] = qttRDVDAC3(d, ori_phi_tt, "TTrounding", bool_energy, bool_save);
    time_TTrounding(i, :) = time;
    round_time_TTrounding(i, :) = sum(round_time);
    energy_TTrounding{i} = energy;
      try
        errors_TTrounding(i, :) = norm(phi - phi_pcg) / norm(phi_pcg);
      catch
        errors_TTrounding(i, j) = -1;
      end
%     errors_TTrounding(i, :) = norm(phi - phi_pcg) / norm(phi_pcg);
  catch
    time_TTrounding(i, :) = inf;
    round_time_TTrounding(i, :) = inf;
  end
  for j = 1:S
    disp(['S = ', num2str(j), ' is being tested.']);
    try
      [round_time, time, phi, energy] = qttRDVDAC3(d, ori_phi_tt, "randorth", bool_energy, bool_save);
      time_randorth(i, j) = time;
      round_time_randorth(i, j) = sum(round_time);
      energy_randorth{i, j} = energy;
      try
        errors_randorth(i, j) = norm(phi - phi_pcg) / norm(phi_pcg);
      catch
        errors_randorth(i, j) = -1;
      end
    catch
      time_randorth(i, j) = inf;
      round_time_randorth(i, j) = inf;
    end
    try
      [round_time, time, phi, energy] = qttRDVDAC3(d, ori_phi_tt, "orthrand", bool_energy, bool_save);
      time_orthrand(i, j) = time;
      round_time_orthrand(i, j) = sum(round_time);
      energy_orthrand{i, j} = energy;
      try
        errors_orthrand(i, j) = norm(phi - phi_pcg) / norm(phi_pcg);
      catch
        errors_orthrand(i, j) = -1;
      end
    catch
      time_orthrand(i, j) = inf;
      round_time_orthrand(i, j) = inf;
    end
    try
      [round_time, time, phi, energy] = qttRDVDAC3(d, ori_phi_tt, "twosided", bool_energy, bool_save);
      time_twosided(i, j) = time;
      round_time_twosided(i, j) = sum(round_time);
      energy_twosided{i, j} = energy;
      try
        errors_twosided(i, j) = norm(phi - phi_pcg) / norm(phi_pcg);
      catch
        errors_twosided(i, j) = -1;
      end
    catch
      time_twosided(i, j) = inf;
      round_time_twosided(i, j) = inf;
    end
    try
      [round_time, time, phi, energy] = qttRDVDAC3(d, ori_phi_tt, "HaTT1", bool_energy, bool_save);
      time_HaTT1(i, j) = time;
      round_time_HaTT1(i, j) = sum(round_time);
      energy_HaTT1{i, j} = energy;
      try
        errors_HaTT1(i, j) = norm(phi - phi_pcg) / norm(phi_pcg);
      catch
        errors_HaTT1(i, j) = -1;
      end
    catch
      time_HaTT1(i, j) = inf;
      round_time_HaTT1(i, j) = inf;
    end
    try
      [round_time, time, phi, energy] = qttRDVDAC3(d, ori_phi_tt, "HaTT2", bool_energy, bool_save);
      time_HaTT2(i, j) = time;
      round_time_HaTT2(i, j) = sum(round_time);
      energy_HaTT2{i, j} = energy;
      try
        errors_HaTT2(i, j) = norm(phi - phi_pcg) / norm(phi_pcg);
      catch
        errors_HaTT2(i, j) = -1;
      end
    catch
      time_HaTT2(i, j) = inf;
      round_time_HaTT2(i, j) = inf;
    end
  end
end
save("qttRDVDAC3.mat", "-regexp", "d_set", "^errors_", "^time_", "^round_time_", "^energy_");