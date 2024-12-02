clear; clc
time_ori = zeros(1, 9);
time_randorth = time_ori;
time_hadamard = time_ori;
D = 4 : 10;
R = 2 .^ (D - 3);
for k = 1 : 7
  d = D(k);
  n = 2^d;
  A = fn(n, n, R(k));
  B = fn(n, n, R(k));
  sz = 2 * ones(1, 2 * d);
  A = reshape(A, sz);
  B = reshape(B, sz);
  TTA = tt_tensor(A);
  TTB = tt_tensor(B);
  l = R(k);
  ori = tic;
  try
    TT = TTA .* TTB;
    round(TT, l);
    time_ori(k) = toc(ori);
  catch
    time_ori(k) = inf;
  end
  randorth = tic;
  try
    TT = TTA .* TTB;
    round_randorth(TT, l);
    time_randorth(k) = toc(randorth);
  catch
    time_randorth(k) = inf;
  end
  hadamard = tic;
  try
    hadamard_round_randorth(TTA, TTB, l);
    time_hadamard(k) = toc(hadamard);
  catch
    time_hadamard(k) = inf;
  end
end
plot(1:7, time_ori, 1:7, time_randorth, 1:7, time_hadamard);
legend('ori', 'randorth', 'hadamard')
