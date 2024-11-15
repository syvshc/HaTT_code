function [round_time, time, phi] = qttRDVDAC3(d, ori_phi, method)
  t0 = cputime;
%   time = 0;
  I = [1, 0; 0, 1];
  J = [0, 1; 0, 0];
  
  n = 2 ^ d;
  hx = 2 * pi / n;
  T = 0.1;
  dt = 0.01;
  gamma = 1;
  beta = 1;
  ep = 0.1;

  sz = 2 * ones(1, 3 * d);
  ori_phi = reshape(ori_phi, sz);
  phi = tt_tensor(ori_phi);

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
  
  clearvars coresB_right coresB_left coresB_mid cores_right cores_mid cores_left coresB_left coresB_mid coresB_right cores_laplace;

  % I = tt_ones(2, d);
  round_time = zeros(T/dt, 1);
  for t = 1:(T/dt)
      aa = reshape(full(phi), n, n, n);
      mesh(aa(:, :, n/2))
      if t == 1
          phi_half = phi;
      else
          phi_half = .5 * (3 * phi - phi1);
      end
      phi_half = round(phi_half, 1e-5);
      test_rank = phi_half.r;
      test_rank(2:end-1) = test_rank(2:end-1) + 3;
      switch method
        case "TTrounding"
          tic;
          E_half = round(phi_half .* phi_half, 1e-5);
          E_half = round(E_half .* phi_half, 1e-5); 
          round_time(t) = toc;
          enddisp = "> TTrounding process ended";
        case "randorth"
          tic;
          E_half = round_randorth(phi_half .* phi_half, test_rank);
          E_half = round_randorth(E_half .* phi_half, test_rank);
          round_time(t) = toc;
          enddisp = "> randorth process ended";
        case "orthrand"
          tic;
          E_half = round_orthrand(phi_half .* phi_half, test_rank);
          E_half = round_orthrand(E_half .* phi_half, test_rank);
          round_time(t) = toc;
          enddisp = "> orthrand process ended";
        case "twosided"
          tic;
          E_half = round_twosided(phi_half .* phi_half, test_rank);
          E_half = round_twosided(E_half .* phi_half, test_rank);
          round_time(t) = toc;
          enddisp = "> twosided process ended";
        case "HaTT1"
          tic;
          E_half = HaTT1(phi_half, phi_half, test_rank);
          E_half = HaTT1(E_half, phi_half, test_rank);
          round_time(t) = toc;
          enddisp = "> HaTT1 process ended";
        case "HaTT2"
          tic;
          E_half = HaTT2(phi_half, phi_half, test_rank);
          E_half = HaTT2(E_half, phi_half, test_rank);
          round_time(t) = toc;
          enddisp = "> HaTT2 process ended";
      end
      E_half = 1 / ep ^ 2 * (E_half - (1 + beta) * phi_half);
      b = B * phi - dt * E_half;
      b = round(b, 1e-5);
      phi1 = phi;
      phi = dmrg_solve2(A, b, 1e-5, 'verb', 0);
      phi = round(phi, 1e-5);
      % if mod(t, 10) == 0
      % %     phi_mat = full(phi);
      % %     phi_mat = reshape(phi_mat, [n, n]);
      % %     contour(xx, xx, phi_mat);
      % %     title(["t=", t]);
      % %     pause;
      %   E_half.r
      % end
      % aa = reshape(full(phi), n, n, n);
      
      % mesh(aa(:, :, n/2))
  end
  time = cputime - t0;
  disp(enddisp)
end