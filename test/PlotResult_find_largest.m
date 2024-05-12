function f = PlotResult_find_largest(file)
    load(file, "err*", "D_set", "errors*", "time*");
    %% Plot results
    
    f = figure(1);
    % f.Position(1:2) = [0,1050];
    % f.Position(3:4) = [1050, 700];
    % f.Position = [-1698.5,-341,1252,1026];
    % Post-process errors
    
    [err_TTrounding, neg_TTrounding, pos_TTrounding] = computeError(errors_TTrounding);
    [err_randorth, neg_randorth, pos_randorth] = computeError(errors_randorth);
    [err_orthrand, neg_orthrand, pos_orthrand] = computeError(errors_orthrand);
    [err_twosided, neg_twosided, pos_twosided] = computeError(errors_twosided);
    [err_HRO, neg_HRO, pos_HRO] = computeError(errors_HRO);
    [err_HRO_no_svd, neg_HRO_no_svd, pos_HRO_no_svd] = computeError(errors_HRO_no_svd);

    subplot(2, 2, 1)
    
    errorbar(D_set, err_TTrounding, neg_TTrounding, pos_TTrounding,     'bo','markersize',10,'linewidth',2)
    ax = gca;
    ax.YScale = 'log';
    
    hold on
    
    errorbar(D_set, err_randorth, neg_randorth, pos_randorth, 'g^','markersize',10,'linewidth',2)
    errorbar(D_set, err_orthrand, neg_orthrand, pos_orthrand, 'm+','markersize',10,'linewidth',2)
    errorbar(D_set, err_twosided, neg_twosided, pos_twosided, 'c*','markersize',10,'linewidth',2)
    errorbar(D_set, err_HRO, neg_HRO, pos_HRO,     'ks','markersize',10,'linewidth',2)
    errorbar(D_set, err_HRO_no_svd, neg_HRO_no_svd, pos_HRO_no_svd, 'rx','markersize',10,'linewidth',2)
    
    hold off
    
    title('(a)')
    xlabel('TT ranks', 'FontSize', 18)
    ylabel('Relative Error', 'FontSize', 18)
    legend('TT-Rounding', 'RandOrth', 'OrthRand', 'TwoSided', 'HRO-1', 'HRO-2', 'Location','southeast')
    set(gca,'FontSize',14,"FontName", "Times New Roman")
    axis square;
    
    % Post-process timings
    
    [speedup_TTrounding, neg_TTrounding, pos_TTrounding] = computeSpeedup(time_TTrounding, time_TTrounding);
    [speedup_randorth, neg_randorth, pos_randorth] = computeSpeedup(time_randorth, time_TTrounding);
    [speedup_orthrand, neg_orthrand, pos_orthrand] = computeSpeedup(time_orthrand, time_TTrounding);
    [speedup_twosided, neg_twosided, pos_twosided] = computeSpeedup(time_twosided, time_TTrounding);
    [speedup_HRO, neg_HRO, pos_HRO] = computeSpeedup(time_HRO, time_TTrounding);
    [speedup_HRO_no_svd, neg_HRO_no_svd, pos_HRO_no_svd] = computeSpeedup(time_HRO_no_svd, time_TTrounding);
    
    subplot(2, 2, 2)
    
    errorbar(D_set, speedup_TTrounding,       neg_TTrounding,       pos_TTrounding,       'bo-','markersize',10,'linewidth',2);
    
    hold on
    
    errorbar(D_set, speedup_randorth,  neg_randorth,  pos_randorth,  'g^-','markersize',10,'linewidth',2);
    errorbar(D_set, speedup_orthrand,  neg_orthrand,  pos_orthrand,  'm+-','markersize',10,'linewidth',2);
    errorbar(D_set, speedup_twosided, neg_twosided, pos_twosided, 'c*-','markersize',10,'linewidth',2);
    errorbar(D_set, speedup_HRO,       neg_HRO,       pos_HRO,       'ks-','markersize',10,'linewidth',2);
    errorbar(D_set, speedup_HRO_no_svd, neg_HRO_no_svd, pos_HRO_no_svd, 'rx-','markersize',10,'linewidth',2);
    
    hold off
    title("(b)")
    xlabel('TT ranks', 'FontSize', 18)
    ylabel('Speedup', 'FontSize', 18)
    % legend('TTrounding', 'HRO', 'HRO-no-SVD', 'randorth', 'orthrand')
    set(gca,'FontSize',14,"FontName", "Times New Roman")
    axis square;

    % Post-process timings

    [times_TTrounding, neg_TTrounding, pos_TTrounding] = computeTime(time_TTrounding);
    [times_randorth, neg_randorth, pos_randorth] = computeTime(time_randorth);
    [times_orthrand, neg_orthrand, pos_orthrand] = computeTime(time_orthrand);
    [times_twosided, neg_twosided, pos_twosided] = computeTime(time_twosided);
    [times_HRO, neg_HRO, pos_HRO] = computeTime(time_HRO);
    [times_HRO_no_svd, neg_HRO_no_svd, pos_HRO_no_svd] = computeTime(time_HRO_no_svd);

    subplot(2, 2, 3)

    errorbar(D_set, times_TTrounding,       neg_TTrounding,       pos_TTrounding,       'bo-','markersize',10,'linewidth',2);

    hold on

    errorbar(D_set, times_randorth,  neg_randorth,  pos_randorth,  'g^-','markersize',10,'linewidth',2);
    errorbar(D_set, times_orthrand,  neg_orthrand,  pos_orthrand,  'm+-','markersize',10,'linewidth',2);
    errorbar(D_set, times_twosided, neg_twosided, pos_twosided, 'c*-','markersize',10,'linewidth',2);
    errorbar(D_set, times_HRO,       neg_HRO,       pos_HRO,       'ks-','markersize',10,'linewidth',2);
    errorbar(D_set, times_HRO_no_svd, neg_HRO_no_svd, pos_HRO_no_svd, 'rx-','markersize',10,'linewidth',2);
    
    hold off
    title("(c)")
    xlabel('TT ranks', 'FontSize', 18)
    ylabel('Time (s)', 'FontSize', 18)
    % legend('TTrounding', 'HRO', 'HRO-no-SVD', 'randorth', 'orthrand')
    set(gca,'FontSize',14,"FontName", "Times New Roman", 'YScale','log')
    axis square;

    % Post-process partical contraction timings
    [times_HBilinearForm, neg_HBilinearForm, pos_HBilinearForm] = computeTime(time_HBilinearForm);
    [times_HBilinearForm_no_svd, neg_HBilinearForm_no_svd, pos_HBilinearForm_no_svd] = computeTime(time_HBilinearForm_no_svd);
    [times_dot, neg_dot, pos_dot] = computeTime(time_dot);

    subplot(2, 2, 4)

    hold on

    errorbar(D_set, times_dot,  neg_dot,  pos_dot,  'g^-','markersize',10,'linewidth',2);
    errorbar(D_set, times_HBilinearForm,       neg_HBilinearForm,       pos_HBilinearForm,       'ks-','markersize',10,'linewidth',2);
    errorbar(D_set, times_HBilinearForm_no_svd, neg_HBilinearForm_no_svd, pos_HBilinearForm_no_svd, 'rx-','markersize',10,'linewidth',2);
    
    hold off
    title("(d)")
    xlabel('TT ranks', 'FontSize', 18)
    ylabel('Time (s)', 'FontSize', 18)
    legend('dot', 'HBilinearForm-1', 'HBilinearForm-2')
    set(gca,'FontSize',14,"FontName", "Times New Roman")
    box on
    axis square;
end