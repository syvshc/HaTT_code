function [X] = round_twosided(Y, varargin)

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
    L = TTrandn(n, l);
    % reset l to a d+1 length vector
    l = L.r;
    rho = [1; floor(2*l(2:d)); 1];
    R = TTrandn(n, rho);
    rho = R.r;
    r_x = l;
    % init ps and core of X
    ps_x = cumsum([1; n(1 : end) .* l(1 : end - 1) .* l(2 : end)]);
    cr_x = zeros(ps_x(end) - 1, 1);
    % generate the sketch phase
    % tic;
    pos = 1;

    % WL = PartialContractionsRL(L, Y);
    % WR = PartialContractionsLR(Y, R);
    WL = PartialContractionsLR(L, Y);
    % WL{k} is a matrix of size r(k) * l(k)
    WR = PartialContractionsRL(Y, R);
    % WR{k} is a matrix of size rho(k) * r(k)
    ps_wr = WR.ps; ps_wl = WL.ps;
    cr_wr = WR.core; cr_wl = WL.core;
    core0 = cr_y(ps_y(1) : ps_y(2) - 1);
    WL0 = cr_wl(ps_wl(1) : ps_wl(2) - 1);
    WR0 = cr_wr(ps_wr(1) : ps_wr(2) - 1);
    core0 = reshape(core0, n(1), r(2));
    WL0 = reshape(WL0, l(2), r(2));
    WR0 = reshape(WR0, r(2), rho(2));
    [U, S, V] = svd(WL0 * WR0, 'econ');
    S = diag(diag(S).^(-0.5));
    L = WR0 * V * S;
    R = S * U' * WL0;

    core0 = core0 * L;
    num0 = numel(core0);
    cr_x(pos : pos + num0 - 1) = core0(:);
    pos = pos + num0;

    for k = 2 : d - 1
      core0 = cr_y(ps_y(k) : ps_y(k + 1) - 1);
      WL0 = cr_wl(ps_wl(k) : ps_wl(k + 1) - 1);
      WR0 = cr_wr(ps_wr(k) : ps_wr(k + 1) - 1);
      core0 = reshape(core0, r(k) * n(k), r(k + 1));
      WL0 = reshape(WL0, l(k + 1), r(k + 1));
      WR0 = reshape(WR0, r(k + 1), rho(k + 1));
      [U, S, V] = svd(WL0 * WR0, 'econ');
      S = diag(diag(S).^(-0.5));
      L = WR0 * V * S;
      core0 = R * v2h(core0 * L, n(k));
      r_x(k) = size(core0, 1);
      num0 = numel(core0);
      cr_x(pos : pos + num0 - 1) = core0(:);
      pos = pos + num0;

      R = S * U' * WL0;
    end
    core0 = cr_y(ps_y(end - 1) : ps_y(end) - 1);
    core0 = reshape(core0, r(end - 1), n(end));
    core0 = R * core0;
    r_x(end - 1) = size(core0, 1);
    num0 = numel(core0);
    cr_x(pos : pos + num0 - 1) = core0(:);
    pos = pos + num0;
    X.core = cr_x(1 : pos - 1);
    X.ps = cumsum([1; n(1 : end) .* r_x(1 : end - 1) .* r_x(2 : end)]);
    X.r = r_x;
    X.n = n;
    X.d = d;
end