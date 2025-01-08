function W = HPCRL1(X, Y, Z)
% PARTIALCONTRACTIONSRL  Compute the sequence of nested contractions between two TT-tensors from right to left, the first TT-tensor is generated by Hadamard product, i.e.do partial contractions of X .* Y and Z
%
%  W = PARTIALCONTRACTIONSRL(X, Y, Z) returns a vector W where reshape(W(ps(i):ps(i + 1) - 1), [X.d(i) * Y.d(i), Z.d(i)] ) is the contraction between cores (i+1, ..., N) of X .* Y and Z.
% \pkp mean partial kronecker product of tensors
  % test if X and Y are TT-tensors
  if (~isa(X, 'tt_tensor') || ~isa(Y, 'tt_tensor') || ~isa(Z, 'tt_tensor'))
    error('X, Y and Z should be TT-tensors');
  end
  % test if X and Y can be Hadamard product together and Z can be contracted with the result
  d_x = X.d; d_y = Y.d; d_z = Z.d;
  n_x = X.n; n_y = Y.n; n_z = Z.n;
  if (any([d_x, d_y, d_z] ~= [d_x, d_x, d_x]) || any([n_x, n_y, n_z] ~= [n_x, n_x, n_x], 'all'))
    error('X, Y and Z should have the same dimensions and sizes');
  else
    d = d_x; 
    n = n_x;
  end
  r_x = X.r; r_y = Y.r; r_z = Z.r;
  ps_x = X.ps; ps_y = Y.ps; ps_z = Z.ps;
  cr_x = X.core; cr_y = Y.core; cr_z = Z.core;
  % initialize the output W
  % W.core is a vector contains all the elements of the W_i
  % W.ps is a vector contains the position of W_i, the same meaning as tt.ps
  W.core = zeros(1, sum(r_x(2 : end - 1) .* r_y(2 : end - 1) .* r_z(2 : end - 1)));
  W.ps = cumsum([1; r_x(2 : end - 1) .* r_y(2 : end - 1) .* r_z(2 : end - 1)]);
  % compute W_d = H(X_d) * H(Y_d)'
  core_x = cr_x(ps_x(end - 1) : ps_x(end) - 1);
  core_x = reshape(core_x, r_x(end - 1), n(end));
  core_y = cr_y(ps_y(end - 1) : ps_y(end) - 1);
  core_y = reshape(core_y, r_y(end - 1), n(end));
  core_z = cr_z(ps_z(end - 1) : ps_z(end) - 1);
  core_z = reshape(core_z, r_z(end - 1), n(end));
  W1 = zeros(r_x(end - 1) * r_y(end - 1), n(end));
  for i = 1:n(end)
    W1(:, i) = kron(core_x(:, i), core_y(:, i));
  end
  W1 = W1 * core_z';
  % put W_d into W.core
  W.core(W.ps(end - 1) : W.ps(end) - 1) = W1(:);
  for k = d - 1 : -1 : 2
    W0 = W1;
    % W0 is a r_y(k+1)r_x(k+1) * r_z(k+1) tensor
    % compute W_(k - 1) = H((X \pkp Y)_(k:n))H(Z_(k:n))':
    core_x = cr_x(ps_x(k) : ps_x(k + 1) - 1);
    core_x = reshape(core_x, r_x(k), n(k), r_x(k + 1));
    core_y = cr_y(ps_y(k) : ps_y(k + 1) - 1);
    core_y = reshape(core_y, r_y(k), n(k), r_y(k + 1));
    core_z = cr_z(ps_z(k) : ps_z(k + 1) - 1);
    core_z = reshape(core_z, r_z(k), n(k), r_z(k + 1));
    W_R = permute(core_z, [1, 3, 2]);
    W_R = reshape(W_R, r_z(k), r_z(k + 1) * n(k));
    W_L = zeros(r_y(k) * r_x(k), n(k) * r_z(k + 1));
    for i = 1 : n(k)
      X1 = reshape(core_x(:, i, :), r_x(k), r_x(k + 1));
      Y1 = reshape(core_y(:, i, :), r_y(k), r_y(k + 1));
      % core_y(:, i, :) is a r_y(k) * 1 * r_y(k + 1) tensor, we reshape it to a matrix
      for j = 1 : r_z(k + 1)
        w = reshape(W0(:, j), r_y(k + 1), r_x(k + 1));
        W_L(:, (i - 1) * r_z(k + 1) + j) = reshape(Y1 * w * X1', [], 1);
      end
    end
    W1 = W_L * W_R';
    % W1 = v2h(core_x * W0, n(k)) * core_y';
    W.core(W.ps(k - 1) : W.ps(k) - 1) = W1(:);
  end
end