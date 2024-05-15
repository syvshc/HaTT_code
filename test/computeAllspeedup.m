function [SU1, SU2] = computeAllspeedup(file)
  load(file, "time*");
  % len = size(time_HaTT, 1);
  % SU1 = zeros(len, 3);
  % SU2 = zeros(len, 3);
  [speedup_randorth, ~, ~] = computeSpeedup(time_HaTT, time_randorth);
  [speedup_orthrand, ~, ~] = computeSpeedup(time_HaTT, time_orthrand);
  [speedup_twosided, ~, ~] = computeSpeedup(time_HaTT, time_twosided);
  SU1 = [speedup_randorth, speedup_orthrand, speedup_twosided];

  [speedup_randorth, ~, ~] = computeSpeedup(time_HaTT_no_svd, time_randorth);
  [speedup_orthrand, ~, ~] = computeSpeedup(time_HaTT_no_svd, time_orthrand);
  [speedup_twosided, ~, ~] = computeSpeedup(time_HaTT_no_svd, time_twosided);
  SU2 = [speedup_randorth, speedup_orthrand, speedup_twosided];
end