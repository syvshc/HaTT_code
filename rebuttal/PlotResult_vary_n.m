load("sincos_n_15.mat", "err*", "test_ranks", "errors*", "time*");
%% Plot results

f = figure('Name', 'vary_n');
% f.Position(1:2) = [0,1050];
% f.Position(3:4) = [1050, 700];
f.Position = [1,49,1536,741.6];
% Post-process errors

[err_TTrounding, neg_TTrounding, pos_TTrounding] = computeError(errors_TTrounding);
[err_randorth, neg_randorth, pos_randorth] = computeError(errors_randorth);
[err_orthrand, neg_orthrand, pos_orthrand] = computeError(errors_orthrand);
[err_twosided, neg_twosided, pos_twosided] = computeError(errors_twosided);
[err_HaTT1, neg_HaTT1, pos_HaTT1] = computeError(errors_HaTT1);
[err_HaTT2, neg_HaTT2, pos_HaTT2] = computeError(errors_HaTT2);

subplot(1, 4, 1)

errorbar(test_ranks, err_TTrounding, neg_TTrounding, pos_TTrounding,     'o', 'Color', '#3570b6','markersize',6,'linewidth',1.5)
% ax = gca;
% ax.YScale = 'log';

hold on

errorbar(test_ranks, err_randorth, neg_randorth, pos_randorth, '^', 'Color', '#23a6ba','markersize',6,'linewidth',1.5)
errorbar(test_ranks, err_orthrand, neg_orthrand, pos_orthrand, '+', 'Color', '#b482ba','markersize',6,'linewidth',1.5)
errorbar(test_ranks, err_twosided, neg_twosided, pos_twosided, '*', 'Color', '#c48f00','markersize',6,'linewidth',1.5)
errorbar(test_ranks, err_HaTT1, neg_HaTT1, pos_HaTT1,     'ks','markersize',6,'linewidth',1.5)
errorbar(test_ranks, err_HaTT2, neg_HaTT2, pos_HaTT2, 'x', 'Color', '#c45c30','markersize',6,'linewidth',1.5)

hold off

title('(1a) $n = 15$','interpreter','latex')
set(gca,'FontSize',14,"FontName", "Times New Roman")
xlabel('TT ranks', 'FontSize', 18)
ylabel('Relative Error (\times 10^{-5})', 'FontSize', 18)
xlim([0, 64])
legend('TT-Rounding', 'RandOrth', 'OrthRand', 'TwoSided', 'HaTT-1', 'HaTT-2')
axis square;


% Post-process timings

[speedup_TTrounding, ~, ~] = computeSpeedup(time_TTrounding, time_TTrounding);
[speedup_randorth, neg_randorth, pos_randorth] = computeSpeedup(time_randorth, time_TTrounding);
[speedup_orthrand, neg_orthrand, pos_orthrand] = computeSpeedup(time_orthrand, time_TTrounding);
[speedup_twosided, neg_twosided, pos_twosided] = computeSpeedup(time_twosided, time_TTrounding);
[speedup_HaTT1, neg_HaTT1, pos_HaTT1] = computeSpeedup(time_HaTT1, time_TTrounding);
[speedup_HaTT2, neg_HaTT2, pos_HaTT2] = computeSpeedup(time_HaTT2, time_TTrounding);

subplot(1, 4, 2)

errorbar(test_ranks, speedup_TTrounding,       zeros(length(test_ranks), 1),       zeros(length(test_ranks), 1),       'o-', 'Color', '#3570b6','markersize',6,'linewidth',1.5);
ax = gca;
ax.YScale = 'log';
ax.XScale = 'log';
hold on

errorbar(test_ranks, speedup_randorth,  neg_randorth,  pos_randorth,  '^-', 'Color', '#23a6ba','markersize',6,'linewidth',1.5);
errorbar(test_ranks, speedup_orthrand,  neg_orthrand,  pos_orthrand,  '+-', 'Color', '#b482ba','markersize',6,'linewidth',1.5);
errorbar(test_ranks, speedup_twosided, neg_twosided, pos_twosided, '*-', 'Color', '#c48f00','markersize',6,'linewidth',1.5);
errorbar(test_ranks, speedup_HaTT1,       neg_HaTT1,       pos_HaTT1,       'ks-','markersize',6,'linewidth',1.5);
errorbar(test_ranks, speedup_HaTT2, neg_HaTT2, pos_HaTT2, 'x-', 'Color', '#c45c30','markersize',6,'linewidth',1.5);
% x = 60:10:110;
% y4 = x.^4/(60^4/speedup_HaTT1(1));
% y3 = x.^3/(60^3/speedup_HaTT1(1));
% y2 = x.^2/(60^2/speedup_HaTT1(1));
% y1 = x/(60/speedup_HaTT1(1));
% % f4 = loglog(x, y4, '--', 'linewidth',1.5);
% f3 = loglog(x, y3, '--', 'linewidth',1.5);
% f2 = loglog(x, y2, '--');
% f1 = loglog(x, y1, '--', 'linewidth', 1.5);
% % legend([f1, f3], ["linear", "cubic"])
% legend([f1, f2, f3], ["linear", "quadratic", "cubic"])
hold off
title("(1b) $n=15$",'interpreter','latex')
set(gca,'FontSize',14,"FontName", "Times New Roman")
xlabel('TT ranks', 'FontSize', 18)
ylabel('Speedup', 'FontSize', 18)
xlim([0, 64])
% legend('TTrounding', 'HaTT', 'HaTT-no-SVD', 'randorth', 'orthrand')

axis square;

load("sincos_n_15.mat", "err*", "test_ranks", "errors*", "time*");

%% Plot results

[err_TTrounding, neg_TTrounding, pos_TTrounding] = computeError(errors_TTrounding);
[err_randorth, neg_randorth, pos_randorth] = computeError(errors_randorth);
[err_orthrand, neg_orthrand, pos_orthrand] = computeError(errors_orthrand);
[err_twosided, neg_twosided, pos_twosided] = computeError(errors_twosided);
[err_HaTT1, neg_HaTT1, pos_HaTT1] = computeError(errors_HaTT1);
[err_HaTT2, neg_HaTT2, pos_HaTT2] = computeError(errors_HaTT2);

subplot(1, 4, 3)

errorbar(test_ranks, err_TTrounding, neg_TTrounding, pos_TTrounding,     'o', 'Color', '#3570b6','markersize',6,'linewidth',1.5)
% ax = gca;
% ax.YScale = 'log';

hold on

errorbar(test_ranks, err_randorth, neg_randorth, pos_randorth, '^', 'Color', '#23a6ba','markersize',6,'linewidth',1.5)
errorbar(test_ranks, err_orthrand, neg_orthrand, pos_orthrand, '+', 'Color', '#b482ba','markersize',6,'linewidth',1.5)
errorbar(test_ranks, err_twosided, neg_twosided, pos_twosided, '*', 'Color', '#c48f00','markersize',6,'linewidth',1.5)
errorbar(test_ranks, err_HaTT1, neg_HaTT1, pos_HaTT1,     'ks','markersize',6,'linewidth',1.5)
errorbar(test_ranks, err_HaTT2, neg_HaTT2, pos_HaTT2, 'x', 'Color', '#c45c30','markersize',6,'linewidth',1.5)

hold off

title('(2a) $n = 20$','interpreter','latex')
set(gca,'FontSize',14,"FontName", "Times New Roman")
xlabel('TT ranks', 'FontSize', 18)
ylabel('Relative Error (\times 10^{-5})', 'FontSize', 18)
xlim([0, 64])
legend('TT-Rounding', 'RandOrth', 'OrthRand', 'TwoSided', 'HaTT-1', 'HaTT-2')
axis square;


% Post-process timings

[speedup_TTrounding, ~, ~] = computeSpeedup(time_TTrounding, time_TTrounding);
[speedup_randorth, neg_randorth, pos_randorth] = computeSpeedup(time_randorth, time_TTrounding);
[speedup_orthrand, neg_orthrand, pos_orthrand] = computeSpeedup(time_orthrand, time_TTrounding);
[speedup_twosided, neg_twosided, pos_twosided] = computeSpeedup(time_twosided, time_TTrounding);
[speedup_HaTT1, neg_HaTT1, pos_HaTT1] = computeSpeedup(time_HaTT1, time_TTrounding);
[speedup_HaTT2, neg_HaTT2, pos_HaTT2] = computeSpeedup(time_HaTT2, time_TTrounding);

subplot(1, 4, 4)

errorbar(test_ranks, speedup_TTrounding,       zeros(length(test_ranks), 1),       zeros(length(test_ranks), 1),       'o-', 'Color', '#3570b6','markersize',6,'linewidth',1.5);
ax = gca;
ax.YScale = 'log';
ax.XScale = 'log';
hold on

errorbar(test_ranks, speedup_randorth,  neg_randorth,  pos_randorth,  '^-', 'Color', '#23a6ba','markersize',6,'linewidth',1.5);
errorbar(test_ranks, speedup_orthrand,  neg_orthrand,  pos_orthrand,  '+-', 'Color', '#b482ba','markersize',6,'linewidth',1.5);
errorbar(test_ranks, speedup_twosided, neg_twosided, pos_twosided, '*-', 'Color', '#c48f00','markersize',6,'linewidth',1.5);
errorbar(test_ranks, speedup_HaTT1,       neg_HaTT1,       pos_HaTT1,       'ks-','markersize',6,'linewidth',1.5);
errorbar(test_ranks, speedup_HaTT2, neg_HaTT2, pos_HaTT2, 'x-', 'Color', '#c45c30','markersize',6,'linewidth',1.5);
% x = 60:10:110;
% y4 = x.^4/(60^4/speedup_HaTT1(1));
% y3 = x.^3/(60^3/speedup_HaTT1(1));
% y2 = x.^2/(60^2/speedup_HaTT1(1));
% y1 = x/(60/speedup_HaTT1(1));
% % f4 = loglog(x, y4, '--', 'linewidth',1.5);
% f3 = loglog(x, y3, '--', 'linewidth',1.5);
% f2 = loglog(x, y2, '--');
% f1 = loglog(x, y1, '--', 'linewidth', 1.5);
% % legend([f1, f3], ["linear", "cubic"])
% legend([f1, f2, f3], ["linear", "quadratic", "cubic"])
hold off
title("(2b) $n=20$",'interpreter','latex')
set(gca,'FontSize',14,"FontName", "Times New Roman")
xlabel('TT ranks', 'FontSize', 18)
ylabel('Speedup', 'FontSize', 18)
xlim([0, 64])
% legend('TTrounding', 'HaTT', 'HaTT-no-SVD', 'randorth', 'orthrand')

axis square;