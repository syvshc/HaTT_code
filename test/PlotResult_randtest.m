function f = PlotResult_randtest(varargin)
  if nargin == 1
    file = varargin{1};
    name = 'randtest';
  else
    file = varargin{1};
    name = varargin{2};
  end
  load(file, "err*", "test_ranks", "errors*", "time*");
  %% Plot results
  
  f = figure('Name', name);
  % f.Position(1:2) = [0,1050];
  % f.Position(3:4) = [1050, 700];
  f.Position = [1,49,1440,781.5];
  % Post-process errors
  
  [err_HaTT, neg_HaTT, pos_HaTT] = computeError(errors_HaTT);
  [err_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd] = computeError(errors_HaTT_no_svd);

  subplot(1, 2, 1)


  hold on

  ax = gca;
  ax.YScale = 'log';
  
  errorbar(test_ranks, err_HaTT, neg_HaTT, pos_HaTT,     'ks','markersize',6,'linewidth',1.5)
  errorbar(test_ranks, err_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd, 'x', 'Color', '#c45c30','markersize',6,'linewidth',1.5)
  
  hold off
  
  title('(a)')
  xlabel('Target ranks', 'FontSize', 24)
  ylabel('Relative Error', 'FontSize', 24)
  legend('HaTT-1', 'HaTT-2')
  set(gca,'FontSize',14,"FontName", "Times New Roman")
  box on
  axis square;

  % Post-process timings

  [times_HaTT, neg_HaTT, pos_HaTT] = computeTime(time_HaTT);
  [times_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd] = computeTime(time_HaTT_no_svd);

  subplot(1, 2, 2)

  hold on

  errorbar(test_ranks, times_HaTT,       neg_HaTT,       pos_HaTT,       'ks-','markersize',6,'linewidth',1.5);
  errorbar(test_ranks, times_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd, 'x-', 'Color', '#c45c30','markersize',6,'linewidth',1.5);
  
  hold off
  title("(c)")
  xlabel('Target ranks', 'FontSize', 24)
  ylabel('Time (s)', 'FontSize', 24)
  % legend('TTrounding', 'HaTT', 'HaTT-no-SVD', 'randorth', 'orthrand')
  set(gca,'FontSize',14,"FontName", "Times New Roman")
  box on
  axis square;
end