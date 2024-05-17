% sincos = PlotResult_sincos('sincos.mat');
% compare12 = PlotResult_compare12('compare12.mat');
% vary_r = PlotResult_vary_r('randtest_vary_r.mat');
qing = PlotResult_find_largest('qing.mat', 'qing');
alpine = PlotResult_find_largest('alpine.mat', 'alpine');

% exportgraphics(sincos, '../fig/sincos.pdf')
% exportgraphics(compare12, '../fig/compare12.pdf')
% exportgraphics(vary_r, '../fig/randtest_vary_r.pdf')
exportgraphics(qing, '../fig/qing.pdf')
exportgraphics(alpine, '../fig/alpine.pdf')