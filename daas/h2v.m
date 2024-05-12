function core = h2v(core, n)
% %V2H change a core from horizontal to vertical format
  core = reshape(core, size(core, 1) * n, []);
end