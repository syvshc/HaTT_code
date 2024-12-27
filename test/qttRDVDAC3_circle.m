function [round_time, time, phi, energy] = qttRDVDAC3_circle(d, ori_phi, method, bool_energy, bool_save)
  t0 = tic;
  energy = 0; 
%   time = 0;
  I = [1, 0; 0, 1];
  J = [0, 1; 0, 0];
  
  n = 2 ^ d;
  hx = 2 / n;
  T = 1;
  dt = 0.01;
ep = 3.5 / 2^(d - 1);
gamma = 6e-5;
% ep = 0.1;
beta = 2;
  
  sz = 2 * ones(1, 3 * d);
  if class(ori_phi) == "double"
    ori_phi = reshape(ori_phi, sz);
  end

  phi = tt_tensor(ori_phi);

  if bool_energy
    r = zeros(1, floor(T/dt) + 1);
    E0 = r;
    % II = reshape(ones(1, n ^ 3), sz);
    % II = tt_tensor(II);
    II = tt_tensor;
    II.core = ones(1, 2 * 3 * d);
    II.r = ones(1, 3 * d + 1);
    II.n = 2 * ones(1, 3 * d);
    II.d = length(II.n);
    II.ps = cumsum([1;II.n.*II.r(1:end - 1).*II.r(2:end)]);
  end

  cores_laplace = cell(d, 1);
  for key = 1:d
    if key == 1
      cores_laplace{key} = zeros(2, 2, 3);
      cores_laplace{key}(:, :, 1) = 1/3 * (1 + 0.5 * dt * gamma * (6 / hx ^ 2 + beta / ep ^ 2)) * I - 0.5 * dt * gamma / hx ^ 2 * J - 0.5 * dt * gamma / hx ^ 2 * J';
      cores_laplace{key}(:, :, 2) = -0.5 * dt * gamma / hx ^ 2 * J;
      cores_laplace{key}(:, :, 3) = -0.5 * dt * gamma / hx ^ 2 * J';
    elseif key == d
      cores_laplace{key} = zeros(2, 2, 3);
      cores_laplace{key}(:, :, 1) = I;
      cores_laplace{key}(:, :, 2) = J';
      cores_laplace{key}(:, :, 3) = J;
    else
      cores_laplace{key} = zeros(2, 2, 3, 3);
      cores_laplace{key}(:, :, 1, 1) = I;
      cores_laplace{key}(:, :, 2, 2) = J;
      cores_laplace{key}(:, :, 3, 3) = J';
      cores_laplace{key}(:, :, 2, 1) = J';
      cores_laplace{key}(:, :, 3, 1) = J;
    end
  end
  cores_laplace{d-1}(:, :, 3, 2) = J';
  cores_laplace{d-1}(:, :, 2, 3) = J;
    
  cores_right = cell(d, 1);
  for key = 1 : d
    if key == 1
      cores_right{key} = zeros(2, 2, 2, 4);
      cores_right{key}(:, :, 1, 1:3) = cores_laplace{1};
      cores_right{key}(:, :, 2, 4) = I;
    elseif key == d
      cores_right{key} = zeros(2, 2, 4);
      cores_right{key}(:, :, 1:3) = cores_laplace{d};
      cores_right{key}(:, :, 4) = I;
    else
      cores_right{key} = zeros(2, 2, 4, 4);
      cores_right{key}(:, :, 1:3, 1:3) = cores_laplace{key};
      cores_right{key}(:, :, 4, 4) = I;
    end
  end

  cores_left = cell(d, 1);
  for key = 1 : d
    if key == 1
      cores_left{key} = zeros(2, 2, 4);
      cores_left{key}(:, :, 1) = I;
      cores_left{key}(:, :, 2:4) = cores_laplace{1};
    elseif key == d
      cores_left{key} = zeros(2, 2, 4, 2);
      cores_left{key}(:, :, 1, 1) = I;
      cores_left{key}(:, :, 2:4, 2) = cores_laplace{d};
    else
      cores_left{key} = zeros(2, 2, 4, 4);
      cores_left{key}(:, :, 1, 1) = I;
      cores_left{key}(:, :, 2:4, 2:4) = cores_laplace{key};
    end
  end

  cores_mid = cell(d, 1);
  for key = 1 : d
    if key == 1
      cores_mid{key} = zeros(2, 2, 2, 5);
      cores_mid{key}(:, :, 1, 1) = I;
      cores_mid{key}(:, :, 1, 2:4) = cores_laplace{1};
      cores_mid{key}(:, :, 2, 5) = I;
    elseif key == d
      cores_mid{key} = zeros(2, 2, 5, 2);
      cores_mid{key}(:, :, 1, 1) = I;
      cores_mid{key}(:, :, 2:4, 2) = cores_laplace{d};
      cores_mid{key}(:, :, 5, 2) = I;
    else
      cores_mid{key} = zeros(2, 2, 5, 5);
      cores_mid{key}(:, :, 1, 1) = I;
      cores_mid{key}(:, :, 2:4, 2:4) = cores_laplace{key};
      cores_mid{key}(:, :, 5, 5) = I;
    end
  end

  
  coresB_right = cores_right;
  coresB_left = cores_left;
  coresB_mid = cores_mid;

  coresB_right{1}(:, :, 1, 1) = 1/3 * (1 - 0.5 * dt * gamma * (6 / hx ^ 2 + beta / ep ^ 2)) * I + 0.5 * dt * gamma / hx ^ 2 * J + 0.5 * dt * gamma / hx ^ 2 * J';
  coresB_right{1}(:, :, 1, 2) = 0.5 * dt * gamma / hx ^ 2 * J;
  coresB_right{1}(:, :, 1, 3) = 0.5 * dt * gamma / hx ^ 2 * J';

  coresB_left{1}(:, :, 2) = 1/3 * (1 - 0.5 * dt * gamma * (6 / hx ^ 2 + beta / ep ^ 2)) * I + 0.5 * dt * gamma / hx ^ 2 * J + 0.5 * dt * gamma / hx ^ 2 * J';
  coresB_left{1}(:, :, 3) = 0.5 * dt * gamma / hx ^ 2 * J;
  coresB_left{1}(:, :, 4) = 0.5 * dt * gamma / hx ^ 2 * J';

  coresB_mid{1}(:, :, 1, 2) = 1/3 * (1 - 0.5 * dt * gamma * (6 / hx ^ 2 + beta / ep ^ 2)) * I + 0.5 * dt * gamma / hx ^ 2 * J + 0.5 * dt * gamma / hx ^ 2 * J';
  coresB_mid{1}(:, :, 1, 3) = 0.5 * dt * gamma / hx ^ 2 * J;
  coresB_mid{1}(:, :, 1, 4) = 0.5 * dt * gamma / hx ^ 2 * J';
  
  A = tt_matrix([cores_left; cores_mid; cores_right]);
  B  = tt_matrix([coresB_left; coresB_mid; coresB_right]);
  
  if bool_energy
    coresL_left = cores_left;
    coresL_mid = cores_mid;
    coresL_right = cores_right;

    coresL_right{1}(:, :, 1, 1) = 1/3 * (6 / hx ^ 2 + beta / ep ^ 2) * I -1 / hx ^ 2 * J -1 / hx ^ 2 * J';
    coresL_right{1}(:, :, 1, 2) = -1 / hx ^ 2 * J;
    coresL_right{1}(:, :, 1, 3) = -1 / hx ^ 2 * J';
  
    coresL_left{1}(:, :, 2) = 1/3 * (6 / hx ^ 2 + beta / ep ^ 2) * I -1 / hx ^ 2 * J -1 / hx ^ 2 * J';
    coresL_left{1}(:, :, 3) = -1 / hx ^ 2 * J;
    coresL_left{1}(:, :, 4) = -1 / hx ^ 2 * J';
  
    coresL_mid{1}(:, :, 1, 2) = 1/3 * (6 / hx ^ 2 + beta / ep ^ 2) * I -1 / hx ^ 2 * J -1 / hx ^ 2 * J';
    coresL_mid{1}(:, :, 1, 3) = -1 / hx ^ 2 * J;
    coresL_mid{1}(:, :, 1, 4) = -1 / hx ^ 2 * J';

    L  = tt_matrix([coresL_left; coresL_mid; coresL_right]);
    r_tmp = round(phi .* phi, 1e-5);
    r_tmp = (1 + beta) * II - r_tmp;
    r(1) = hx * .25 / ep ^ 2 * dot(r_tmp, r_tmp);
    E0(1) = .5 * gamma * hx * dot(phi, L * phi);
  end

  clearvars coresB_right coresB_left coresB_mid cores_right cores_mid cores_left coresB_left coresB_mid coresB_right cores_laplace coresL_left coresL_mid coresL_right r_tmp II;  

  % I = tt_ones(2, d);
  round_time = zeros(T/dt, 1);
  for t = 1:(T/dt)
    if bool_save
      eval(['phi_hatt2_', num2str(t), '= reshape(full(phi), n, n, n);']);
      if isfile("phi_hatt2.mat") == 1
        eval(['save(''phi_hatt2.mat'', ''phi_hatt2_', num2str(t), ''', ''-append'')']);
      else
        eval(['save(''phi_hatt2.mat'', ''phi_hatt2_', num2str(t), ''')']);
      end
    end
      if t == 1
          phi_half = phi;
      else
          phi_half = .5 * (3 * phi - phi1);
      end
      phi_half = round(phi_half, 1e-5);
      test_rank = phi_half.r;
      test_rank(2:end-1) = test_rank(2:end-1);
      switch method
        case "TTrounding"
          tic;
          E_half = round(phi_half .* phi_half, 1e-5);
          E_half = round(E_half .* phi_half, 1e-5); 
          round_time(t) = toc;
          enddisp = ['> d = ', num2str(d), ': TTrounding process ended'];
        case "randorth"
          tic;
          E_half = round_randorth(phi_half .* phi_half, test_rank);
          E_half = round_randorth(E_half .* phi_half, test_rank);
          round_time(t) = toc;
          enddisp = ['> d = ', num2str(d), ': randorth process ended'];
        case "orthrand"
          tic;
          E_half = round_orthrand(phi_half .* phi_half, test_rank);
          E_half = round_orthrand(E_half .* phi_half, test_rank);
          round_time(t) = toc;
          enddisp = ['> d = ', num2str(d), ': orthrand process ended'];
        case "twosided"
          tic;
          E_half = round_twosided(phi_half .* phi_half, test_rank);
          E_half = round_twosided(E_half .* phi_half, test_rank);
          round_time(t) = toc;
          enddisp = ['> d = ', num2str(d), ': twosided process ended'];
        case "HaTT1"
          tic;
          E_half = HaTT1(phi_half, phi_half, test_rank);
          E_half = HaTT1(E_half, phi_half, test_rank);
          round_time(t) = toc;
          enddisp = ['> d = ', num2str(d), ': HaTT1 process ended'];
        case "HaTT2"
          tic;
          E_half = HaTT2(phi_half, phi_half, test_rank);
          E_half = HaTT2(E_half, phi_half, test_rank);
          round_time(t) = toc;
          enddisp = ['> d = ', num2str(d), ': HaTT2 process ended'];
      end
      E_half = 1 / ep ^ 2 * (E_half - (1 + beta) * phi_half);
      b = B * phi - dt * E_half;
      b = round(b, 1e-5);
      phi1 = phi;
      phi = dmrg_solve2(A, b, 1e-5, 'verb', 0);
      phi = round(phi, 1e-5);
      if mod(t, 10) == 0
          phi_mat = full(phi);
          phi_mat = reshape(phi_mat, [n, n, n]);
          mesh(phi_mat(:, :, n/2));
          title(["t=", t]);
          % pause;
      %   E_half.r
      end
      % aa = reshape(full(phi), n, n, n);
      
      % mesh(aa(:, :, n/2))
      
    if bool_energy
      r(t + 1) = r(t) + hx * dot(E_half, phi - phi1);
      E0(t + 1) = .5 * gamma * hx * dot(phi, L * phi);
    end
  end
  if bool_energy  
    energy = E0 + r;
  end
  time = toc(t0);
  disp(enddisp)
end