function W = PartialContractionsLR(X, Y)
  % PARTIALCONTRACTIONSRL  Compute the sequence of nested contractions between two TT-tensors from right to left.
  %
  %  W = PARTIALCONTRACTIONSRL(X, Y) returns a vector W where reshape(W(ps(i):ps(i + 1) - 1), [X.d(i), Y.d(i)] ) is the contraction between cores (i+1, ..., N) of X and Y.
  %
    % test if X and Y are TT-tensors
    if (~isa(X, 'tt_tensor') || ~isa(Y, 'tt_tensor'))
      error('X and Y should be TT-tensors');
    end
    % test if X and Y can be contracted together
    d_x = X.d; d_y = Y.d;
    n_x = X.n; n_y = Y.n;
    if (any(d_x ~= d_y) || any(n_x ~= n_y))
      error('X and Y should have the same dimensions and sizes');
    else
      d = d_x; 
      n = n_x;
    end
    r_x = X.r; r_y = Y.r;
    ps_x = X.ps; ps_y = Y.ps;
    cr_x = X.core; cr_y = Y.core;
    % initialize the output W
    % W.core is a vector contains all the elements of the W_i
    % W.ps is a vector contains the position of W_i, the same meaning as tt.ps
    W.core = zeros(1, sum(r_x(2 : end - 1).*r_y(2 : end - 1)));
    W.ps = cumsum([1; r_x(2 : end - 1) .* r_y(2 : end - 1)]);
    % compute W_1 = V(X_1)' * V(Y_1)
    core_x = cr_x(ps_x(1) : ps_x(2) - 1);
    core_x = reshape(core_x, n(1), r_x(2));
    core_y = cr_y(ps_y(1) : ps_y(2) - 1);
    core_y = reshape(core_y, n(2), r_y(2));
    W1 = core_x' * core_y;
    % put W_1 into W.core
    W.core(W.ps(1) : W.ps(2) - 1) = W1(:);
    for k = 2 : d - 1
      W0 = W1;
      % compute W_(k) = V(X_(1:k))'V(Y_(1:k)):
      % H(Z_k) = W_k*H(X_k)*W_(k-1)
      % W_(k) = V(Z_k)'V(Y_k)
      core_x = cr_x(ps_x(k) : ps_x(k + 1) - 1);
      core_x = reshape(core_x, r_x(k) * n(k), r_x(k + 1));
      core_y = cr_y(ps_y(k) : ps_y(k + 1) - 1);
      core_y = reshape(core_y, r_y(k),  n(k) * r_y(k + 1));
      W1 = core_x' * h2v(W0 * core_y, n(k));
      W.core(W.ps(k) : W.ps(k + 1) - 1) = W1(:);
    end
  end