function f = PlotResult_compare12(varargin)
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
  f.Position = [1,49,1536,741.6];
  % Post-process errors
  
  [err_HaTT1, neg_HaTT1, pos_HaTT1] = computeError(errors_HaTT1);
  [err_HaTT2, neg_HaTT2, pos_HaTT2] = computeError(errors_HaTT2);

  subplot(1, 2, 1)


  hold on

  ax = gca;
  ax.YScale = 'log';
  
  errorbar(test_ranks, err_HaTT1, neg_HaTT1, pos_HaTT1,     'ks','markersize',6,'linewidth',1.5)
  errorbar(test_ranks, err_HaTT2, neg_HaTT2, pos_HaTT2, 'x', 'Color', '#c45c30','markersize',6,'linewidth',1.5)
  
  hold off
  
  title('(a)')
  set(gca,'FontSize',14.625,"FontName", "Times New Roman")
  xlabel('Target ranks', 'FontSize', 21.375)
  ylabel('Relative Error (\times 10^{-15})', 'FontSize', 21.375)
  xlim([25, 125])
  legend('HaTT-1', 'HaTT-2', 'Location','southeast')
  box on
  axis square;

  % Post-process timings

  [times_HaTT1, neg_HaTT1, pos_HaTT1] = computeTime(time_HaTT1);
  [times_HaTT2, neg_HaTT2, pos_HaTT2] = computeTime(time_HaTT2);

  subplot(1, 2, 2)

  hold on

  errorbar(test_ranks, times_HaTT1,       neg_HaTT1,       pos_HaTT1,       'ks-','markersize',6,'linewidth',1.5);
  errorbar(test_ranks, times_HaTT2, neg_HaTT2, pos_HaTT2, 'x-', 'Color', '#c45c30','markersize',6,'linewidth',1.5);
  
  hold off
  title("(b)")
  set(gca,'FontSize',14.625,"FontName", "Times New Roman")
  xlabel('Target ranks', 'FontSize', 21.375)
  ylabel('Time (s)', 'FontSize', 21.375)
  xlim([25, 125])
  % legend('TTrounding', 'HaTT', 'HaTT-no-SVD', 'randorth', 'orthrand')
  box on
  axis square;
end