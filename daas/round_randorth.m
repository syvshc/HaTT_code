function [varargout] = round_randorth(Y, varargin)
% ROUND_RANDORTH  Round a TT-tensor to a TT-tensor with randomization then orthogonalization.
% Y is a TT-tensor, 
% If the target TT-rank is not given, we'll use the size of each unfolding of Y as the default target TT-rank.

  % init X
  X = tt_tensor;
  n = Y.n; d = Y.d; r = Y.r;
  cr_y = Y.core; ps_y = Y.ps;
  % if l is given, use it, otherwise, generate l
  % we generate l as the max rank of Y's each unfolding
  if nargin == 2
    l = varargin{1};
  elseif nargin == 1
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
  time_Mcore0 = 0;
  % generate the sketch phase
  PCRL = tic;
  W = PartialContractionsRL(Y, R);
  time_PCRL = toc(PCRL);
  % time_tmp = toc;
  % disp(["time of PartialContractionsRL for Y.r = ", num2str(Y.r')])
  % disp(["R.r = ", num2str(l')])
  % disp(["time = ", num2str(time_tmp)])
  cr_w = W.core; ps_w = W.ps;
  % init ps and core of X
  ps_x = cumsum([1; n(1 : end) .* l(1 : end - 1) .* l(2 : end)]);
  cr_x = zeros(ps_x(end) - 1, 1);
  % set V(Y_1) to core0
  core0 = cr_y(ps_y(1) : ps_y(2) - 1);
  core0 = reshape(core0, [r(1) * n(1), r(2)]);
  pos = 1;
  for i = 1 : d - 1
    % set W_i to core1
    % save core0 to Z
    core1 = cr_w(ps_w(i) : ps_w(i + 1) - 1);
    core1 = reshape(core1, [r(i + 1), l(i + 1)]);
    Z = core0;
    % compute V(X_i)W_i and econ-QR on it
    % set Q to core0 i.e. V(X_i)
    [core0, ~] = qr(Z * core1, 0);
    % M = V(X_i)'*Z
    M = core0' * Z;
    % renew r_x(i + 1) if l(i + 1) > r(i + 1)
    r_x(i + 1) = size(core0, 2);
    num0 = numel(core0);
    % put V(X_i) to cr_x
    cr_x(pos : pos + num0 - 1) = core0(:);
    pos = pos + num0;
    % reest H(Y_(i+1)) to core0
    % and renew H(X_(i+1)) = M*H(Y_(i+1))
    core0 = cr_y(ps_y(i + 1) : ps_y(i + 2) - 1);
    core0 = reshape(core0, [r(i + 1), n(i + 1) * r(i + 2)]);
    Mcore0 = tic;
    core0 = M * core0;
    time_Mcore0 = time_Mcore0 + toc(Mcore0);
    core0 = h2v(core0, n(i + 1));
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
      varargout{2} = time_PCRL;
      varargout{3} = time_Mcore0;
        
    case 2 
      varargout{1} = X;
      varargout{2} = time_PCRL;
    case 1
      varargout{1} = X;  
      
  end
end