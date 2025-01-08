function [time, phi] = RDVDAC2(d)
  tic
  % n = 2 ^ d;
  % hx = 2 * pi / n;
  % dt = 0.01;
  % beta = 1;
  % T = 1;
  % ep = 0.1;
  % gamma = 0.5;
  % xx = ((1:n)' - 0.5) * hx;
% %   phi = sin(xx);
%   phi = ones(n, 1);
% %   phi = 0.2 * (phi * phi');
%   phi = phi * phi';
  % phi = reshape(ori_phi, [n ^ 2, 1]);
  % circle
  n = 2 ^ d;
  hx = 256 / n;
  R0 = 100;
  ep = 0.0078;
  % gamma = 1;
  % ep = 0.1;
  gamma = 6e-5;
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

  T = 0.5;
  dt = 0.01;
  beta = 0;
  phi = reshape(ori_phi, [n ^ 2, 1]);

  % generate 1-D sparse A and B and a identity matrix I
  i_AB = [1:n, 1:(n - 1), 2:n, n, 1];
  j_AB = [1:n, 2:n, 1:(n - 1),  1, n];
  v_A = [0.5 * (1 + 0.5 * dt * gamma * (4 / hx ^ 2 + beta / ep ^ 2)) * ones(1, n), -0.5 * dt * gamma / hx ^ 2 * ones(1, 2 * n)];
  v_B = [0.5 * (1 - 0.5 * dt * gamma * (4 / hx ^ 2 + beta / ep ^ 2)) * ones(1, n), 0.5 * dt * gamma / hx ^ 2 * ones(1, 2 * n)];
  A = sparse(i_AB, j_AB, v_A);
  B = sparse(i_AB, j_AB, v_B);
  
  i_I = 1:n;
  v_I = ones(1, n);
  I = sparse(i_I, i_I, v_I);

  % generate 2-D sparse A and B
  A = kron(A, I) + kron(I, A);
  B = kron(B, I) + kron(I, B);
  L = ichol(A);
  for t = 1 : T/dt
    if t == 1
        phi_half = phi;
    else
        phi_half = 0.5 * (3 * phi - phi1);
    end
    b = B * phi - dt * gamma * phi_half .* (phi_half .* phi_half - 1 - beta) / ep^2;
    phi1 = phi; % 存储 phi_{n-1} 的值
    [phi, flag] = pcg(A, b);
    % phi = A\b;
    % phiphi = reshape(phi, [n, n]); mesh(phiphi);
%     phiphi = reshape(phi, [n, n]);
%     mesh(xx, xx, phiphi)
%     pause
      % phi_tt = reshape(phi, sz);
      % phi_tt = tt_tensor(phi_tt);
      % phi_tt = round(phi_tt, 1e-5);
      % r = phi_tt.r;
    %   mesh(reshape(phi, [n, n]))
  end
  % phi = reshape(phi, [n, n]);
  % mesh(xx, xx, phi);
  time = toc;
end