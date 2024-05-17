%% Custom function to post-process errors

function [error, neg, pos] = computeError(errors)
  errors = sort(errors, 2);
  errors = errors(:,2:end-1);
  error = squeeze(mean(errors, 2));
  neg = squeeze(std(errors,[],2));
  pos = squeeze(std(errors,[],2));
end