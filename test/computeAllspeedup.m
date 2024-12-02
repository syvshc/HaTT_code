function [SU1, SU2] = computeAllspeedup(file)
% function [TIME1] = computeAllspeedup(file)
  load(file, "time*");
  % len = size(time_HaTT1, 1);
  % SU1 = zeros(len, 3);
  % SU2 = zeros(len, 3);
  [speedup_randorth, ~, ~] = computeSpeedup(time_HaTT1, time_randorth);
  [speedup_orthrand, ~, ~] = computeSpeedup(time_HaTT1, time_orthrand);
  [speedup_twosided, ~, ~] = computeSpeedup(time_HaTT1, time_twosided);
  % [time_HaTT1, ~, ~] = computeTime(time_HaTT1);
  % [time_HaTT2, ~, ~] = computeTime(time_HaTT2);
  % [time_randorth, ~, ~] = computeTime(time_randorth);
  % [time_orthrand, ~, ~] = computeTime(time_orthrand);
  % [time_twosided, ~, ~] = computeTime(time_twosided);
  % [time_TTrounding, ~, ~] = computeTime(time_TTrounding);
  % TIME1 = [time_HaTT1, time_HaTT2, time_randorth, time_orthrand, time_twosided, time_TTrounding];
  SU1 = [speedup_randorth, speedup_orthrand, speedup_twosided];

  [speedup_randorth, ~, ~] = computeSpeedup(time_HaTT2, time_randorth);
  [speedup_orthrand, ~, ~] = computeSpeedup(time_HaTT2, time_orthrand);
  [speedup_twosided, ~, ~] = computeSpeedup(time_HaTT2, time_twosided);
  SU2 = [speedup_randorth, speedup_orthrand, speedup_twosided];
end