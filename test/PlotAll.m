sincos = PlotResult_sincos('sincos.mat');
randtest = PlotResult_randtest('randtest.mat');
vary_r = PlotResult_vary_r('randtest_vary_r.mat');
qing = PlotResult_find_largest('qing.mat', 'qing');
alpine = PlotResult_find_largest('alpine.mat', 'alpine');

exportgraphics(sincos, '../fig/sincos.pdf')
exportgraphics(randtest, '../fig/randtest.pdf')
exportgraphics(vary_r, '../fig/randtest_vary_r.pdf')
exportgraphics(qing, '../fig/qing.pdf')
exportgraphics(alpine, '../fig/alpine.pdf')