load("phi_hatt2.mat")
load("phi_pcg.mat")
f = figure("Name","DrawAC");
hx = 2*pi/2^7;
xx = hx * (1:2^7);
t_set = 2:2:10;
label_set = {'(a)', '(b)', '(c)', '(d)', '(e)'};
contour_set = {'[-.3 -.2 -.1 .1 .2 .3]',...
                '[-.8 -.6 -.4 -.2 .2 .4 .6 .8]',...
                '[-.8 -.6 -.4 -.2 .2 .4 .6 .8]',...
                '[-.8 -.6 -.4 -.2 .2 .4 .6 .8]',...
                '[-.8 -.6 -.4 -.2 .2 .4 .6 .8]'};
colorbar_set = {'-.3:.1:.3', '-.8:.2:.8', '-.8:.2:.8', '-.8:.2:.8', '-.8:.2:.8'};

tiledlayout(2, 3);
for i = 1:5
  nexttile;
  % subplot(2, 3, i);
  eval(['contourf(xx, xx, phi_pcg_', num2str(t_set(i)), '(:, :, 32), ', contour_set{i}, ', ''k-'', "LineWidth", 1);'])
  hold on;
  eval(['contour(xx, xx, phi_hatt2_', num2str(t_set(i)), '(:, :, 32), ', contour_set{i}, ', ''r--'', "LineWidth", 1);'])
  axis square
  set(gca,'xtick',0:pi/2:2*pi);
  set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'});
  set(gca,'ytick',0:pi/2:2*pi);
  set(gca,'YTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'});
  colormap("sky");
  eval(['colorbar(''Ticks'', ', colorbar_set{i}, ')']);
  title([label_set{i}, ' $t = ', num2str(t_set(i)*0.01), '$'], 'Interpreter', 'latex');
  k = plot(NaN, 'k-', "LineWidth", 1.5);
  r = plot(NaN, 'r--', "LineWidth", 1.5);
end
  lgd = legend([k, r], {'PCG', 'HaTT2'}, 'FontSize',12);
  lgd.Layout.Tile = 6;