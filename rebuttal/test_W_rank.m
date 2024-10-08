%% This is a maniscript for testing if W has low-rank structure
clear

n = 20 * ones(7, 1);
d = length(n);
r = [1; 100 * ones(d - 1, 1); 1];

%% Set cores

y = tt_tensor();
y.n = n; y.d = d; y.r = r;
y.ps = cumsum([1; r(1 : end-1) .* n .* r(2 : end)]);
y.core = zeros(y.ps(end)-1, 1);
for k = 1 : d
  core = zeros(r(k), n(k), r(k + 1));
  for i1 = 1 : r(k)
    for i2 = 1 : n(k)
      for i3 = 1 : r(k + 1)
        core(i1, i2, i3) = 1 / (i1 + i2 + i3 - 1);
      end
    end
  end
  y.core(y.ps(k) : y.ps(k + 1) - 1) = core(:);
end
z = y;
%%
% item_num = 60;
% d = 7;
% n = 10;
% len = n^d;
% t1 = 2*pi/len:2*pi/len:2*pi;
% a = 10 * rand(1, item_num) + 0.1;
% b = 10 * rand(1, item_num) + 0.1;
% y1 = 0; z1 = 0;
% for k = 1: item_num
%   y1 = y1 + a(k) * sin(k * t1);
%   z1 = z1 + b(k) * cos(k * t1);
% end
% y = y1; z = z1;

% n = n * ones(1, d);
% y = tt_tensor(reshape(y, n))
% z = tt_tensor(reshape(z, n))
% clear y1 z1 t1 a b;

%%

% ranks to test every 2
test_ranks = [20, 40, 60];

L = length(test_ranks);
% all Wk is of size 100 * l
% S_size = zeros(L, S);
% S_size_trunc = zeros(L, S);
f = figure('Name','testWsv');
f.Position = [46.6,199.4,1404.8,416];
for i = 1 : L
    ell = test_ranks(i);
    R = TTrandn(n, ell);
    l = R.r;
    W_trunc = HPCRL1(y, z, R);
    for k = 1 : d - 1
      W1 = W_trunc.core(W_trunc.ps(k) : W_trunc.ps(k + 1) - 1);
      W1 = reshape(W1, [], l(k + 1));
      if k == 2
          [U, S, V] = svd(W1, 'econ'); 
          subplot(1, 3, i);
          semilogy(diag(S));
          axis square;
          ylim([1e-20, 1])
          title(['target rank = ', num2str(ell)])
      end
    end
end