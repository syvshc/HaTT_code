function z = HBF2(A, varargin)
  % This function computes the inner product of x and A.*y for nargin == 3 or A.*x for nargin == 2;
  % A, X, y are all vectors in tensor train format.
  % In matrix form, we consider matrix B as a diagonal matrix with A on the diagonal
  % so we get a Biliear form of the form x'By or x'Bx
  % We take advantage of the structure of Hadamard product A.*y to accelarate the computation.

  if nargin == 3
    x = varargin{1};
    y = varargin{2};
  elseif nargin == 2
    x = varargin{1};
    y = x;
  end
  if ~isa(x, 'tt_tensor') || ~isa(y, 'tt_tensor') || ~isa(A, 'tt_tensor');
    error("Invalid input arguments.") 
  end
  % test if A and y can be Hadamard product together and x can be contracted with the result
  d_x = x.d; d_y = y.d; d_a = A.d;
  n_x = x.n; n_y = y.n; n_a = A.n;
  if (any([d_x, d_y, d_a] ~= [d_x, d_x, d_x]) || any([n_x, n_y, n_a] ~= [n_x, n_x, n_x], 'all'))
    error('x, y and A should have the same dimensions and sizes');
  else
    d = d_x; 
    n = n_x;
  end
  r_x = x.r; r_y = y.r; r_a = A.r;
  ps_x = x.ps; ps_y = y.ps; ps_a = A.ps;
  cr_x = x.core; cr_y = y.core; cr_a = A.core;
    % compute W_d = H(A_d) * H(y_d)'
    core_x = cr_x(ps_x(end - 1) : ps_x(end) - 1);
    core_x = reshape(core_x, r_x(end - 1), n(end));
    core_y = cr_y(ps_y(end - 1) : ps_y(end) - 1);
    core_y = reshape(core_y, r_y(end - 1), n(end));
    core_a = cr_a(ps_a(end - 1) : ps_a(end) - 1);
    core_a = reshape(core_a, r_a(end - 1), n(end));
    W1 = zeros(r_a(end - 1) * r_y(end - 1), n(end));
    for i = 1:n(end)
      W1(:, i) = kron(core_a(:, i), core_y(:, i));
    end
    W1 = W1 * core_x';
    for k = d - 1 : -1 : 1
      W0 = W1;
      % W1 = zeros(r_a(k) * r_y(k), r_x(k));
      % [U, S, V] = svd(W0, 'econ', 'vector');
      % compute W_(k - 1) = H((A \pkp y)_(k:n))H(x_(k:n))':
      core_x = cr_x(ps_x(k) : ps_x(k + 1) - 1);
      core_x = reshape(core_x, r_x(k), n(k), r_x(k + 1));
      W_R = permute(core_x, [1, 3, 2]);
      W_R = reshape(W_R, r_x(k), r_x(k + 1) * n(k));
      core_y = cr_y(ps_y(k) : ps_y(k + 1) - 1);
      core_y = reshape(core_y, r_y(k), n(k), r_y(k + 1));
      core_a = cr_a(ps_a(k) : ps_a(k + 1) - 1);
      core_a = reshape(core_a, r_a(k), n(k), r_a(k + 1));
      W_L = zeros(r_a(k) * r_y(k), n(k) * r_x(k + 1));
      % W1 = zeros(r_a(k) * r_y(k), r_x(k));
      for i = 1 : n(k)
        % core_y(:, i, :) is a r_y(k) * 1 * r_y(k + 1) tensor, we reshape it to a matrix
        % X1 = reshape(core_x(:, i, :), r_x(k), r_x(k + 1));
        Y1 = reshape(core_y(:, i, :), r_y(k), r_y(k + 1));
        A1 = reshape(core_a(:, i, :), r_a(k), r_a(k + 1));
        for j = 1 : r_x(k + 1)
          w = reshape(W0(:, j), r_y(k + 1), r_a(k + 1));
          W_L(:, (i - 1) * r_x(k + 1) + j) = reshape(Y1 * w * A1', [], 1);
          % W_L(:, j) = reshape(Y1 * w * A1', [], 1);
        end
        % W1 = W1 + W_L * X1';
      end
      W1 = W_L * W_R';
    end
    z = W1;
end



% function z = HBF1(A, varargin)
%   % This function computes the inner product of x and A.*y for nargin == 3 or A.*x for nargin == 2;
%   % A, X, y are all vectors in tensor train format.
%   % In matrix form, we consider matrix B as a diagonal matrix with A on the diagonal
%   % so we get a Biliear form of the form x'By or x'Bx
%   % We take advantage of the structure of Hadamard product A.*y to accelarate the computation.

%   if nargin == 3
%     x = varargin{1};
%     y = varargin{2};
%   elseif nargin == 2
%     x = varargin{1};
%     y = x;
%   end
%   if ~isa(x, 'tt_tensor') || ~isa(y, 'tt_tensor') || ~isa(A, 'tt_tensor');
%     error("Invalid input arguments.") 
%   end
%   % test if A and y can be Hadamard product together and x can be contracted with the result
%   d_x = x.d; d_y = y.d; d_a = A.d;
%   n_x = x.n; n_y = y.n; n_a = A.n;
%   if (any([d_x, d_y, d_a] ~= [d_x, d_x, d_x]) || any([n_x, n_y, n_a] ~= [n_x, n_x, n_x], 'all'))
%     error('x, y and A should have the same dimensions and sizes');
%   else
%     d = d_x; 
%     n = n_x;
%   end
%   r_x = x.r; r_y = y.r; r_a = A.r;
%   ps_x = x.ps; ps_y = y.ps; ps_a = A.ps;
%   cr_x = x.core; cr_y = y.core; cr_a = A.core;
%     % compute W_d = H(A_d) * H(y_d)'
%     core_x = cr_x(ps_x(end - 1) : ps_x(end) - 1);
%     core_x = reshape(core_x, r_x(end - 1), n(end));
%     core_y = cr_y(ps_y(end - 1) : ps_y(end) - 1);
%     core_y = reshape(core_y, r_y(end - 1), n(end));
%     core_a = cr_a(ps_a(end - 1) : ps_a(end) - 1);
%     core_a = reshape(core_a, r_a(end - 1), n(end));
%     W1 = zeros(r_a(end - 1) * r_y(end - 1), n(end));
%     for i = 1:n(end)
%       W1(:, i) = kron(core_a(:, i), core_y(:, i));
%     end
%     W1 = W1 * core_x';
%     for k = d - 1 : -1 : 1
%       W0 = W1;
%       % W1 = zeros(r_a(k) * r_y(k), r_x(k));
%       [U, S, V] = svd(W0, 'econ', 'vector');
%       % compute W_(k - 1) = H((A \pkp y)_(k:n))H(x_(k:n))':
%       core_x = cr_x(ps_x(k) : ps_x(k + 1) - 1);
%       core_x = reshape(core_x, r_x(k), n(k), r_x(k + 1));
%       core_y = cr_y(ps_y(k) : ps_y(k + 1) - 1);
%       core_y = reshape(core_y, r_y(k), n(k), r_y(k + 1));
%       core_a = cr_a(ps_a(k) : ps_a(k + 1) - 1);
%       core_a = reshape(core_a, r_a(k), n(k), r_a(k + 1));
%       W_L = zeros(r_a(k) * r_y(k), n(k) * length(S));
%       W_R = zeros(r_x(k), n(k) * length(S));
%       for i = 1 : n(k)
%         for j = 1 : length(S)
%           u = reshape(U(:, j), r_y(k + 1), r_a(k + 1));
%           v = V(:, j);
%           % core_y(:, i, :) is a r_y(k) * 1 * r_y(k + 1) tensor, we reshape it to a matrix
%           X1 = reshape(core_x(:, i, :), r_x(k), r_x(k + 1));
%           Y1 = reshape(core_y(:, i, :), r_y(k), r_y(k + 1));
%           A1 = reshape(core_a(:, i, :), r_a(k), r_a(k + 1));
%           W_L(:, (i - 1) * length(S) + j) = reshape(Y1 * u * A1', [], 1);
%           W_R(:, (i - 1) * length(S) + j) = X1 * v;
%           % W1 = W1 + S(j) * reshape(core_y(:, i, :) * u * core_a(:, i, :)', [], 1) * (core_x(:, i, :) * v)';
%         end
%       end
%       W1 = W_L * kron(eye(n(k)), diag(S)) * W_R';
%     end
%     z = W1;
% end