function max_energy_error = computeRDVDEnergy(x, y)
  if ~(isequal(size(x), size(y)))
    error('x and y must have the same size')
  end
  % max_energy_error = zeros(1, size(d_set));
  energy_error = (x - y)./y;
  max_energy_error = max(energy_error, [], [1, 3]);
end