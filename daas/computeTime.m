%% Custom function to post-process runtimes

function [times, neg, pos] = computeTime(time)
  time = sort(time, 2);
  time = time(:,2:end-1);
  times = mean(time, 2);
  neg = std(time,[],2);
  pos = std(time,[],2);
end