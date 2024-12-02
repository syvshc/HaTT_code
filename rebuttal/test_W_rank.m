%% This is a maniscript for testing if W has low-rank structure
clear

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
% 
% sz = n * ones(1, d);
% y = tt_tensor(reshape(y, sz))
% z = tt_tensor(reshape(z, sz))
% clear y1 z1 t1 a b;
d = 7;
n = 20 * ones(d, 1);
r = 60;
y = TTrand(n, r);
z = TTrand(n, r);
% ranks to test every 2
test_ranks = [60, 100, 140];

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
          axis square;
          semilogy(diag(S));
          title(['target rank = ', num2str(ell)])
      end
    end
end