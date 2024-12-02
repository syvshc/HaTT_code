clear
d_set = 5:9;

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
  clearvars X Y Z xx;

  [time, phi_pcg, energy] = RDVDAC3(d, ori_phi, bool_energy);
  time_pcg(i, :) = time;
  energy_pcg{i} = energy;
  phi_pcg = reshape(phi_pcg, 2 * ones(1, 3 * d));
  phi_pcg = tt_tensor(phi_pcg);
  try 
    [round_time, time, phi, energy] = qttRDVDAC3(d, ori_phi, "TTrounding", bool_energy);
    time_TTrounding(i, :) = time;
    round_time_TTrounding(i, :) = sum(round_time);
    energy_TTrounding{i} = energy;
    errors_TTrounding(i, :) = norm(phi - phi_pcg) / norm(phi_pcg);
  catch
    time_TTrounding(i, :) = inf;
    round_time_TTrounding(i, :) = inf;
  end
  for j = 1:S
    disp(['S = ', num2str(j), ' is being tested.']);
    try
      [round_time, time, phi, energy] = qttRDVDAC3(d, ori_phi, "randorth", bool_energy);
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
      [round_time, time, phi, energy] = qttRDVDAC3(d, ori_phi, "orthrand", bool_energy);
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
      [round_time, time, phi, energy] = qttRDVDAC3(d, ori_phi, "twosided", bool_energy);
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
      [round_time, time, phi, energy] = qttRDVDAC3(d, ori_phi, "HaTT1", bool_energy);
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
      [round_time, time, phi, energy] = qttRDVDAC3(d, ori_phi, "HaTT2", bool_energy);
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