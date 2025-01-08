function [varargout] = HaTT1(Y, Z, varargin)
% HADMARD_ROUND  Round the Hadamard product of two TT-tensor with RandOrth
% Y and Z are TT-tensors of the same physics modes,
% If the target TT-rank is not given, we'll use the size of each unfolding of Y as the default target TT-rank.

   % init X
  X = tt_tensor;
  if (~isa(Y, 'tt_tensor') || ~isa(Z, 'tt_tensor'))
    error('Y and Z should be TT-tensors');
  end
  % test if Y and Z can be hadamard producted
  d_y = Y.d; d_z = Z.d;
  n_y = Y.n; n_z = Z.n;
  if ((d_y ~= d_z) && (n_y ~= n_z))
    error('Y and Z should have the same dimensions and sizes');
  else
    d = d_y; 
    n = n_y;
  end
  r_y = Y.r; r_z = Z.r;
  cr_y = Y.core; ps_y = Y.ps;
  cr_z = Z.core; ps_z = Z.ps;
  % if l is given, use it, otherwise, generate l
  % we generate l as the max rank of Y's each unfolding

  if nargin == 3
    l = varargin{1};
  elseif nargin == 2
    l = zeros(d + 1, 1);
    l(1) = 1; l(end) = 1;
    for i = 1 : d - 1
      L = prod(n(1 : i));
      R = prod(n(i + 1 : end));
      l(i + 1) = min(L, R);
    end
  else
    error("Invalid number of arguments.");
  end
  % generate a Gaussian tensor train format.
  R = TTrandn(n, l);
  % reset l to a d+1 length vector
  l = R.r;
  r_x = l;
  time_HMcore0_no_svd = 0;
  % generate the sketch phase
  HPCRL_no_svd = tic;
  W = HPCRL1(Y, Z, R);
  time_HPCRL_no_svd = toc(HPCRL_no_svd);
  cr_w = W.core; ps_w = W.ps;
  % init ps and core of X
  ps_x = cumsum([1; n(1 : end) .* l(1 : end - 1) .* l(2 : end)]);
  cr_x = zeros(ps_x(end) - 1, 1);
  % set V(Y_1) to core0
  corey = cr_y(ps_y(1) : ps_y(2) - 1);
  corey = reshape(corey, [r_y(1), n(1), r_y(2)]);
  corey = permute(corey, [1, 3, 2]);
  corez = cr_z(ps_z(1) : ps_z(2) - 1);
  corez = reshape(corez, [r_z(1), n(1), r_z(2)]);
  corez = permute(corez, [1, 3, 2]);
  % core0 = cr_y(ps_y(1) : ps_y(2) - 1);
  core0 = zeros(r_y(1) * r_z(1), r_y(2) * r_z(2), n(1));
  for j = 1 : n(1)
    core0(:, :, j) = kron(corey(:, :, j), corez(:, :, j));
  end
  core0 = permute(core0, [1, 3, 2]);
  core0 = reshape(core0, [r_y(1) * r_z(1) * n(1), r_y(2) * r_z(2)]);
  pos = 1;
  for k = 1 : d - 1
    % set W_k to core1
    % save core0 to ZZ
    core1 = cr_w(ps_w(k) : ps_w(k + 1) - 1);
    core1 = reshape(core1, [r_y(k + 1) * r_z(k + 1), l(k + 1)]);
    ZZ = core0;
    % compute V(X_k)W_k and econ-QR on it
    % set Q to core0 i.e. V(X_k)
    [core0, ~] = qr(ZZ * core1, 0);
    % M = V(X_k)'*ZZ
    M = core0' * ZZ;
    % renew r_x(k + 1) if l(k + 1) > r(k + 1)
    r_x(k + 1) = size(core0, 2);
    num0 = numel(core0);
    % put V(X_k) to cr_x
    cr_x(pos : pos + num0 - 1) = core0(:);
    pos = pos + num0;
    % reest H(Y_(k+1)) to core0
    % and renew H(X_(k+1)) = M*H(Y_(k+1))
    corey = cr_y(ps_y(k + 1) : ps_y(k + 2) - 1);
    corey = reshape(corey, [r_y(k + 1), n(k + 1), r_y(k + 2)]);
    corey = permute(corey, [1, 3, 2]);
    corez = cr_z(ps_z(k + 1) : ps_z(k + 2) - 1);
    corez = reshape(corez, [r_z(k + 1), n(k + 1), r_z(k + 2)]);
    corez = permute(corez, [1, 3, 2]);
    % % core0 = cr_y(ps_y(1) : ps_y(2) - 1);
    % core00 = zeros(r_y(k + 1) * r_z(k + 1), r_y(k + 2) * r_z(k + 2), n(k + 1));
    % for j = 1 : n(k + 1)
    %   core00(:, :, j) = kron(corey(:, :, j), corez(:, :, j));
    % end
    % core00 = permute(core00, [1, 3, 2]);
    % core00 = reshape(core00, [r_y(k + 1) * r_z(k + 1), n(k + 1) * r_y(k + 2) * r_z(k + 2)]);
    % core00 = M * core00;
    % core00 = h2v(core00, n(k + 1));
    HMcore0_no_svd = tic;
    M_row = size(M, 1);
    core0 = zeros(M_row, n(k + 1), r_y(k + 2) * r_z(k + 2));
    for i = 1 : n(k + 1)
      for beta = 1 : M_row
        m = M(beta, :);
        M1 = reshape(m, r_z(k + 1), r_y(k + 1));
        core0(beta, i, :) = reshape(corez(:, :, i)' * M1 * corey(:, :, i), 1, []);
      end
    end
    time_HMcore0_no_svd = time_HMcore0_no_svd + toc(HMcore0_no_svd);
    core0 = reshape(core0, M_row * n(k + 1), r_y(k + 2) * r_z(k + 2));
  end
  % put the last core into cr_x, and truncate useless elements in cr_x
  num0 = numel(core0);
  cr_x(pos : pos + num0 - 1) = core0(:);
  pos = pos + num0;
  X.core = cr_x(1 : pos - 1);
  X.ps = cumsum([1; n(1 : end) .* r_x(1 : end - 1) .* r_x(2 : end)]);
  X.r = r_x;
  X.n = n;
  X.d = d;
  switch nargout
    case 3
      varargout{1} = X;
      varargout{2} = time_HPCRL_no_svd;
      varargout{3} = time_HMcore0_no_svd;
        
    case 2 
      varargout{1} = X;
      varargout{2} = time_HPCRL_no_svd;
    case 1
      varargout{1} = X;  
      
  end
end