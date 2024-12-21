function f = PlotResult_energy(varargin)
  if nargin == 1
      file = varargin{1};
      name = 'energy';
  else
      file = varargin{1};
      name = varargin{2};
  end
  load(file, "d_set", "energy*");
  % load(file)
  %% Plot results
  
  f = figure('Name', name);
  % f.Position(1:2) = [0,1050];
  % f.Position(3:4) = [1050, 700];
  f.Position = [1,49,1440,781.5];
  % Post-process errors
  d = ceil((d_set(end) + d_set(1)) / 2);
  i = 3; j = 4;
  t = 0 : 0.01 : 0.1;
  subplot(1, 2, 1)

  plot(t, energy_mat_pcg(:, i, j), 'p-', 'Color', '#5D31C4','markersize',6,'linewidth',1.5)
  % ax = gca;
  % ax.YScale = 'log';
  
  hold on
  
  plot(t, energy_mat_TTrounding(:, i, j), 'o-', 'Color', '#3570b6','markersize',6,'linewidth',1.5)
  plot(t, energy_mat_randorth(:, i, j), '^-', 'Color', '#23a6ba','markersize',6,'linewidth',1.5)
  plot(t, energy_mat_orthrand(:, i, j), '+-', 'Color', '#b482ba','markersize',6,'linewidth',1.5)
  plot(t, energy_mat_twosided(:, i, j), '*-', 'Color', '#c48f00','markersize',6,'linewidth',1.5)
  plot(t, energy_mat_HaTT1(:, i, j), 'ks-','markersize',6,'linewidth',1.5)
  plot(t, energy_mat_HaTT2(:, i, j), 'x-', 'Color', '#c45c30','markersize',6,'linewidth',1.5)
  
  hold off
  
  title('(a)')
  xlabel('t', 'FontSize', 18)
  ylabel('Free Energy', 'FontSize', 18)
  legend('PCG', 'TT-Rounding', 'RandOrth', 'OrthRand', 'TwoSided', 'HaTT-1', 'HaTT-2')
  set(gca,'FontSize',14,"FontName", "Times New Roman")
  axis square;

  % Post-process timings
  max_energy_error_TTrounding = computeRDVDEnergy(energy_mat_TTrounding, energy_mat_pcg);
  max_energy_error_randorth = computeRDVDEnergy(energy_mat_randorth, energy_mat_pcg);
  max_energy_error_orthrand = computeRDVDEnergy(energy_mat_orthrand, energy_mat_pcg);
  max_energy_error_twosided = computeRDVDEnergy(energy_mat_twosided, energy_mat_pcg);
  max_energy_error_HaTT1 = computeRDVDEnergy(energy_mat_HaTT1, energy_mat_pcg);
  max_energy_error_HaTT2 = computeRDVDEnergy(energy_mat_HaTT2, energy_mat_pcg);

  subplot(1, 2, 2)

  plot(d_set, max_energy_error_TTrounding, 'o', 'Color', '#3570b6','markersize',6,'linewidth',1.5)
  hold on
  plot(d_set, max_energy_error_randorth, '^', 'Color', '#23a6ba','markersize',6,'linewidth',1.5)
  plot(d_set, max_energy_error_orthrand, '+', 'Color', '#b482ba','markersize',6,'linewidth',1.5)
  plot(d_set, max_energy_error_twosided, '*', 'Color', '#c48f00','markersize',6,'linewidth',1.5)
  plot(d_set, max_energy_error_HaTT1, 'ks','markersize',6,'linewidth',1.5)
  plot(d_set, max_energy_error_HaTT2, 'x', 'Color', '#c45c30','markersize',6,'linewidth',1.5)
  
  hold off
  title("(b)")
  xlabel('d', 'FontSize', 18)
  ylabel('RelErr of Free Energy', 'FontSize', 18)
  % legend('TTrounding', 'HaTT', 'HaTT-no-SVD', 'randorth', 'orthrand')
  set(gca,'FontSize',14,"FontName", "Times New Roman", 'YScale','log')
  xlim([4.5, 9.5]);
  xticks(5:9)
  axis square;
end