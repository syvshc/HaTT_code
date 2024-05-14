function f = PlotResult(file)
    load(file, "err*", "test_ranks", "errors*", "time*");
    %% Plot results
    
    f = figure(1);
    % f.Position(1:2) = [0,1050];
    % f.Position(3:4) = [1050, 700];
    % f.Position = [-2014.5,-72,1919.5,592];
    % Post-process errors
    
    [err_TTrounding, neg_TTrounding, pos_TTrounding] = computeError(errors_TTrounding);
    [err_randorth, neg_randorth, pos_randorth] = computeError(errors_randorth);
    [err_orthrand, neg_orthrand, pos_orthrand] = computeError(errors_orthrand);
    [err_HaTT, neg_HaTT, pos_HaTT] = computeError(errors_HaTT);
    [err_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd] = computeError(errors_HaTT_no_svd);
    [err_twosided, neg_twosided, pos_twosided] = computeError(errors_twosided);

    subplot(1, 3, 1)


    errorbar(test_ranks, err_TTrounding, neg_TTrounding, pos_TTrounding,     'bo','markersize',6,'linewidth',1)
    ax = gca;
    ax.YScale = 'log';
    
    hold on
    
    errorbar(test_ranks, err_randorth, neg_randorth, pos_randorth, 'g^','markersize',6,'linewidth',1)
    errorbar(test_ranks, err_orthrand, neg_orthrand, pos_orthrand, 'm+','markersize',6,'linewidth',1)
    errorbar(test_ranks, err_twosided, neg_twosided, pos_twosided, 'c*','markersize',6,'linewidth',1)
    errorbar(test_ranks, err_HaTT, neg_HaTT, pos_HaTT,     'ks','markersize',6,'linewidth',1)
    errorbar(test_ranks, err_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd, 'rx','markersize',6,'linewidth',1)
    
    hold off
    
    % title('(a)')
    xlabel('Target ranks', 'FontSize', 18)
    ylabel('Relative Error', 'FontSize', 18)
    legend('TT-Rounding', 'RandOrth', 'OrthRand', 'TwoSided', 'HaTT-1', 'HaTT-2')
    set(gca,'FontSize',14,"FontName", "Times New Roman")
    axis square;
    
    % Post-process timings
    
    [speedup_HaTT, neg_HaTT, pos_HaTT] = computeSpeedup(time_HaTT, time_TTrounding);
    [speedup_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd] = computeSpeedup(time_HaTT_no_svd, time_TTrounding);
    [speedup_TTrounding, neg_TTrounding, pos_TTrounding] = computeSpeedup(time_TTrounding, time_TTrounding);
    [speedup_randorth, neg_randorth, pos_randorth] = computeSpeedup(time_randorth, time_TTrounding);
    [speedup_orthrand, neg_orthrand, pos_orthrand] = computeSpeedup(time_orthrand, time_TTrounding);
    [speedup_twosided, neg_twosided, pos_twosided] = computeSpeedup(time_twosided, time_TTrounding);
    
    subplot(1, 3, 2)
    
    errorbar(test_ranks, speedup_TTrounding,       neg_TTrounding,       pos_TTrounding,       'bo-','markersize',6,'linewidth',1);
    
    hold on
    
    errorbar(test_ranks, speedup_randorth,  neg_randorth,  pos_randorth,  'g^-','markersize',6,'linewidth',1);
    errorbar(test_ranks, speedup_orthrand,  neg_orthrand,  pos_orthrand,  'm+-','markersize',6,'linewidth',1);
    errorbar(test_ranks, speedup_twosided, neg_twosided, pos_twosided, 'c*-','markersize',6,'linewidth',1);
    errorbar(test_ranks, speedup_HaTT,       neg_HaTT,       pos_HaTT,       'ks-','markersize',6,'linewidth',1);
    errorbar(test_ranks, speedup_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd, 'rx-','markersize',6,'linewidth',1);
    
    hold off
    % title("(b)")
    xlabel('Target ranks', 'FontSize', 18)
    ylabel('Speedup', 'FontSize', 18)
    % legend('TTrounding', 'HaTT', 'HaTT-no-SVD', 'randorth', 'orthrand')
    set(gca,'FontSize',14,"FontName", "Times New Roman")
    axis square;

    % Post-process timings

    [times_TTrounding, neg_TTrounding, pos_TTrounding] = computeTime(time_TTrounding);
    [times_randorth, neg_randorth, pos_randorth] = computeTime(time_randorth);
    [times_orthrand, neg_orthrand, pos_orthrand] = computeTime(time_orthrand);
    [times_twosided, neg_twosided, pos_twosided] = computeTime(time_twosided);
    [times_HaTT, neg_HaTT, pos_HaTT] = computeTime(time_HaTT);
    [times_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd] = computeTime(time_HaTT_no_svd);

    subplot(1, 3, 3)

    errorbar(test_ranks, times_TTrounding,       neg_TTrounding,       pos_TTrounding,       'bo-','markersize',6,'linewidth',1);

    hold on

    errorbar(test_ranks, times_randorth,  neg_randorth,  pos_randorth,  'g^-','markersize',6,'linewidth',1);
    errorbar(test_ranks, times_orthrand,  neg_orthrand,  pos_orthrand,  'm+-','markersize',6,'linewidth',1);
    errorbar(test_ranks, times_twosided, neg_twosided, pos_twosided, 'c*-','markersize',6,'linewidth',1);
    errorbar(test_ranks, times_HaTT,       neg_HaTT,       pos_HaTT,       'ks-','markersize',6,'linewidth',1);
    errorbar(test_ranks, times_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd, 'rx-','markersize',6,'linewidth',1);
    
    hold off
    % title("(c)")
    xlabel('Target ranks', 'FontSize', 18)
    ylabel('Time (s)', 'FontSize', 18)
    % legend('TTrounding', 'HaTT', 'HaTT-no-SVD', 'randorth', 'orthrand')
    set(gca,'FontSize',14,"FontName", "Times New Roman")
    axis square;
end