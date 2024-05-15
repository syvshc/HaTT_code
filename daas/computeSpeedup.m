%% Custom function to post-process runtimes

function [speedup, neg, pos] = computeSpeedup(time, timeref)
  time = sort(time, 2);
  time = time(:,2:end-1);
  timeref = sort(timeref, 2);
  timeref = timeref(:,2:end-1);
  speedup = mean(timeref, 2) ./ mean(time, 2);
  neg = speedup - mean(timeref,2) ./ max(time,[],2);
  pos = mean(timeref,2) ./ min(time,[],2) - speedup;
end