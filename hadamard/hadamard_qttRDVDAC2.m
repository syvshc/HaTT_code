function [round_time, time, phi] = hadamard_qttRDVDAC2(d)
  t0 = cputime;
%   time = 0;
  I = [1, 0; 0, 1];
  J = [0, 1; 0, 0];
  
%   n = 2 ^ d;
%   hx = 2 * pi / n;
%   T = 1;
%   dt = 0.01;
%   gamma = 1;
%   beta = 1;
%   ep = 0.1;
%   xx = hx * ((1:n)' + 0.5);
%   ori_phi = .2 * sin(xx);
%   ori_phi = ori_phi * ori_phi';

  % circle
  n = 2 ^ d;
  hx = 2 / n;
  R0 = 100/128;
  ep = 3.5 / 2^(d - 1);
  gamma = 1;
%   ep = 0.1;
%   gamma = 6e-5;
  xx = ((1:n) - n/2) * hx;
  ori_phi = zeros(n, n);
  for i = 1:n
      for j = 1:n
%           ori_phi(i, j) = -tanh(2 * (sqrt(xx(i) ^ 2 + xx(j) ^ 2) - R0) / ep);
          if xx(i) ^ 2 + xx(j) ^ 2 < R0 ^ 2
              ori_phi(i, j) = 1;
          else
              ori_phi(i, j) = -1;
          end
      end
  end
% ori_phi = (abs(xx)<R0);
% ori_phi = ori_phi' * ori_phi;
% ori_phi(ori_phi == 0) = -1;
%   mesh(ori_phi);
%   hold on;

  T = 1;
  dt = 0.001;
  beta = 2;

  sz = 2 * ones(1, 2 * d);
  ori_phi = reshape(ori_phi, sz);
  phi = tt_tensor(ori_phi);
  
  cores_1 = cell(d, 1);
  for key = 1 : d
    if key == 1
      cores_1{key} = zeros(2, 2, 2, 4);
      cores_1{key}(:, :, 1, 1) = I;
      cores_1{key}(:, :, 2, 2) = 0.5 * (1 + 0.5 * dt * gamma * (4 / hx ^ 2 + beta / ep ^ 2)) * I - 0.5 * dt * gamma / hx ^ 2 * J - 0.5 * dt * gamma / hx ^ 2 * J';
      cores_1{key}(:, :, 2, 3) = -0.5 * dt * gamma / hx ^ 2 * J;
      cores_1{key}(:, :, 2, 4) = -0.5 * dt * gamma / hx ^ 2 * J';
    elseif key == d
      cores_1{key} = zeros(2, 2, 4);
      cores_1{key}(:, :, 1) = I;
      cores_1{key}(:, :, 2) = I;
      cores_1{key}(:, :, 3) = J';
      cores_1{key}(:, :, 4) = J;
    else
      cores_1{key} = zeros(2, 2, 4, 4);
      cores_1{key}(:, :, 1, 1) = I;
      cores_1{key}(:, :, 2, 2) = I;
      cores_1{key}(:, :, 3, 3) = J;
      cores_1{key}(:, :, 4, 4) = J';
      cores_1{key}(:, :, 3, 2) = J';
      cores_1{key}(:, :, 4, 2) = J;
    end
  end
  cores_1{d - 1}(:, :, 4, 3) = J';
  cores_1{d - 1}(:, :, 3, 4) = J;

  cores_2 = cell(d, 1);
  for key = 1 : d
    if key == 1
      cores_2{key} = zeros(2, 2, 4);
      cores_2{key}(:, :, 1) = 0.5 * (1 + 0.5 * dt * gamma * (4 / hx ^ 2 + beta / ep ^ 2)) * I - 0.5 * dt * gamma / hx ^ 2 * J - 0.5 * dt * gamma / hx ^ 2 * J';
      cores_2{key}(:, :, 2) = -0.5 * dt * gamma / hx ^ 2 * J;
      cores_2{key}(:, :, 3) = -0.5 * dt * gamma / hx ^ 2 * J';
      cores_2{key}(:, :, 4) = I;
    elseif key == d
      cores_2{key} = zeros(2, 2, 4, 2);
      cores_2{key}(:, :, 4, 2) = I;
      cores_2{key}(:, :, 1, 1) = I;
      cores_2{key}(:, :, 2, 1) = J';
      cores_2{key}(:, :, 3, 1) = J;
    else
      cores_2{key} = zeros(2, 2, 4, 4);
      cores_2{key}(:, :, end, end) = I;
      cores_2{key}(:, :, 1, 1) = I;
      cores_2{key}(:, :, 2, 2) = J;
      cores_2{key}(:, :, 3, 3) = J';
      cores_2{key}(:, :, 2, 1) = J';
      cores_2{key}(:, :, 3, 1) = J;
    end
  end
  cores_2{d-1}(:, :, 3, 2) = J';
  cores_2{d-1}(:, :, 2, 3) = J;
  
  coresB_1 = cores_1;
  coresB_2 = cores_2;

  coresB_1{1}(:, :, 2, 2) = 0.5 * (1 - 0.5 * dt * gamma * (4 / hx ^ 2 + beta / ep ^ 2)) * I + 0.5 * dt * gamma / hx ^ 2 * J + 0.5 * dt * gamma / hx ^ 2 * J';
  coresB_1{1}(:, :, 2, 3) = 0.5 * dt * gamma / hx ^ 2 * J;
  coresB_1{1}(:, :, 2, 4) = 0.5 * dt * gamma / hx ^ 2 * J';

  coresB_2{1}(:, :, 1) = 0.5 * (1 - 0.5 * dt * gamma * (4 / hx ^ 2 + beta / ep ^ 2)) * I + 0.5 * dt * gamma / hx ^ 2 * J + 0.5 * dt * gamma / hx ^ 2 * J';
  coresB_2{1}(:, :, 2) = 0.5 * dt * gamma / hx ^ 2 * J;
  coresB_2{1}(:, :, 3) = 0.5 * dt * gamma / hx ^ 2 * J';

  
  A = tt_matrix([cores_2; cores_1]);
  B  = tt_matrix([coresB_2; coresB_1]);
  % I = tt_ones(2, d);
  round_time = zeros(T/dt, 3);
  for t = 1:(T/dt)
      if t == 1
          phi_half = phi;
      else
          phi_half = .5 * (3 * phi - phi1);
      end
      r = phi.r;
      % r(2 : end - 1) = r(2 : end - 1);
      tic;
      phi_half = round_randorth(phi_half, r);
      round_time(t, 1) = toc;
      tic;
      E_half = hadamard_round_randorth(phi_half, phi_half, r);
      round_time(t, 2) = toc;
      tic;
      E_half = hadamard_round_randorth(E_half, phi_half, r); 
      round_time(t, 3) = toc;
      E_half = 1 / ep ^ 2 * (E_half - (1 + beta) * phi_half);
      b = B * phi - dt * E_half;
      b = round(b, 1e-5);
      phi1 = phi;
      phi = dmrg_solve2(A, b, 1e-5, 'x0', phi1, 'verb', 0);
      phi = round(phi, 1e-5);
      % if mod(t, 700) == 0
      %     phi_mat = full(phi);
      %     phi_mat = reshape(phi_mat, [n, n]);
      %     contour(xx, xx, phi_mat);
      %     title(["t=", t]);
      %     pause;
      % end
  end
  time = cputime - t0;

end