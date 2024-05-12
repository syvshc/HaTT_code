%% Custom function to post-process errors

function [error, neg, pos] = computeError(errors)
  errors = sort(errors, 2);
  errors = errors(:,2:end-1);
  error = squeeze(median(errors, 2));
  neg = error -  squeeze(min(errors,[],2));
  pos = squeeze(max(errors,[],2)) - error;
end