function [core] = v2h(core, n)
%V2H change a core from vertical to horizontal format
  core = reshape(core, [], n * size(core, 2));
end