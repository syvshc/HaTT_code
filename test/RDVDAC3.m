function [time, phi, energy] = RDVDAC3(d, ori_phi, bool_energy);
  tic
  energy = 0;
  n = 2 ^ d;
  hx = 2 * pi / n;
  dt = 0.01;
  beta = 1;
  T = 0.1;
  ep = 0.1;
  gamma = 1;

  phi = reshape(ori_phi, [n ^ 3, 1]);

  if bool_energy
    r = zeros(1, floor(T/dt) + 1);
    E0 = r;
  end
  % generate 1-D sparse A and B and a identity matrix I
  i_AB = [1:n, 1:(n - 1), 2:n, n, 1];
  j_AB = [1:n, 2:n, 1:(n - 1),  1, n];
  v_A = [1/3 * (1 + 0.5 * dt * gamma * (6 / hx ^ 2 + beta / ep ^ 2)) * ones(1, n), -0.5 * dt * gamma / hx ^ 2 * ones(1, 2 * n)];
  v_B = [1/3 * (1 - 0.5 * dt * gamma * (6 / hx ^ 2 + beta / ep ^ 2)) * ones(1, n), 0.5 * dt * gamma / hx ^ 2 * ones(1, 2 * n)];
  A = sparse(i_AB, j_AB, v_A);
  B = sparse(i_AB, j_AB, v_B);
  
  i_I = 1:n;
  v_I = ones(1, n);
  I = sparse(i_I, i_I, v_I);

  % generate 2-D sparse A and B
  A = kron(kron(A, I), I) + kron(kron(I, A), I) + kron(kron(I, I), A);
  B = kron(kron(B, I), I) + kron(kron(I, B), I) + kron(kron(I, I), B);

  if bool_energy
    v_L = [1/3 * (6 / hx ^ 2 + beta / ep ^ 2) * ones(1, n), -1 / hx ^ 2 * ones(1, 2 * n)];
    LL = sparse(i_AB, j_AB, v_L);
    LL = kron(kron(LL, I), I) + kron(kron(I, LL), I) + kron(kron(I, I), LL);
    r_tmp = (1 + beta) - phi .* phi;
    r(1) = hx * .25 / ep ^ 2 * dot(r_tmp, r_tmp);
    E0(1) = .5 * gamma * hx * dot(phi, LL * phi);
  end

  L = ichol(A);
  for t = 1 : T/dt
    if t == 1 %#ok<ALIGN>
        phi_half = phi;
    else
        phi_half = 0.5 * (3 * phi - phi1);
      end
    E_half = phi_half .* (phi_half .* phi_half - 1 - beta) / ep^2; 
    b = B * phi - dt * E_half;
    phi1 = phi; % 存储 phi_{n-1} 的值
    [phi, ~] = pcg(A, b, 1e-5,100, L, L');

    if bool_energy
      r(t + 1) = r(t) + hx * dot(E_half, phi - phi1);
      E0(t + 1) = .5 * gamma * hx * dot(phi, LL * phi);
    end
  end
  if bool_energy
    energy = E0 + r;
  end
  phi = reshape(phi, [n, n, n]);
  % mesh(xx, xx, phi);
  time = toc;
  disp("PCG process ends");
end