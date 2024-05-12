function tt = TTrand(n, l)
  %TTrand generates a random tensor train of size [n] and rank [l]
  % [n] - size of the tensor
  % [l] - rank of the randomized tensor train. If [l] is a scalar, then all
  %       TT-ranks except the first and the last are equal to [l]. If [l] is a vector of length(n) + 1, then
  %       the TT-ranks are equal to [l],
  
    % Check the input
    % if length(l) == 1, we set all the TT-ranks to l except the first and the last
    if (length(l) == 1 && length(n) > 1)
      l = [1, ones(1, length(n)-1)*l, 1];
      for i = 1 : length(n) - 1
        L = prod(n(1 : i));
        R = prod(n(i + 1 : end));
        l(i + 1) = min([L, R, l(i + 1)]);
      end
    elseif length(l) ~= 1 && (l(1) ~= 1 || l(end) ~= 1)
      error('The first and the last TT-ranks must be equal to 1');
    end
    % set tt's fields, we reshape l to a column vector
    tt = tt_tensor;
    tt.d = length(n);
    tt.n = reshape(n, tt.d, 1);
    tt.r = reshape(l, tt.d + 1, 1);
    tt.ps = cumsum([1; tt.n .* tt.r(1 : end - 1) .* tt.r(2 : end)]);
    tt.core = zeros(tt.ps(end) - 1, 1);
    for i = 1 : tt.d
      rand_tmp = rand(tt.n(i) * tt.r(i) * tt.r(i + 1), 1);
      rand_tmp = (rand_tmp - min(rand_tmp)) / (max(rand_tmp) - min(rand_tmp));
      tt.core(tt.ps(i) : tt.ps(i + 1) - 1) = rand_tmp;
      clear rand_tmp;
    end
  end