function f = PlotResult_find_largest(varargin)
    if nargin == 1
        file = varargin{1};
        name = 'find_largest';
    else
        file = varargin{1};
        name = varargin{2};
    end
    load(file, "err*", "D_set", "errors*", "time*");
    %% Plot results
    
    f = figure('Name', name);
    % f.Position(1:2) = [0,1050];
    % f.Position(3:4) = [1050, 700];
    f.Position = [1,49,1440,781.5];
    % Post-process errors
    
    [err_TTrounding, neg_TTrounding, pos_TTrounding] = computeError(errors_TTrounding);
    [err_randorth, neg_randorth, pos_randorth] = computeError(errors_randorth);
    [err_orthrand, neg_orthrand, pos_orthrand] = computeError(errors_orthrand);
    [err_twosided, neg_twosided, pos_twosided] = computeError(errors_twosided);
    [err_HaTT, neg_HaTT, pos_HaTT] = computeError(errors_HaTT);
    [err_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd] = computeError(errors_HaTT_no_svd);

    subplot(1, 3, 1)
    
    errorbar(D_set, err_TTrounding, neg_TTrounding, pos_TTrounding,     'o', 'Color', '#3570b6','markersize',6,'linewidth',1.5)
    ax = gca;
    ax.YScale = 'log';
    
    hold on
    
    errorbar(D_set, err_randorth, neg_randorth, pos_randorth, '^', 'Color', '#23a6ba','markersize',6,'linewidth',1.5)
    errorbar(D_set, err_orthrand, neg_orthrand, pos_orthrand, '+', 'Color', '#b482ba','markersize',6,'linewidth',1.5)
    errorbar(D_set, err_twosided, neg_twosided, pos_twosided, '*', 'Color', '#c48f00','markersize',6,'linewidth',1.5)
    errorbar(D_set, err_HaTT, neg_HaTT, pos_HaTT,     'ks','markersize',6,'linewidth',1.5)
    errorbar(D_set, err_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd, 'x', 'Color', '#c45c30','markersize',6,'linewidth',1.5)
    
    hold off
    
    % title('Alpine-(a)')
    xlabel('Dimension', 'FontSize', 18)
    ylabel('Relative Error', 'FontSize', 18)
    legend('TT-Rounding', 'RandOrth', 'OrthRand', 'TwoSided', 'HaTT-1', 'HaTT-2', 'Location','southeast')
    set(gca,'FontSize',14,"FontName", "Times New Roman")
    axis square;

    % Post-process timings

    [times_TTrounding, neg_TTrounding, pos_TTrounding] = computeTime(time_TTrounding);
    [times_randorth, neg_randorth, pos_randorth] = computeTime(time_randorth);
    [times_orthrand, neg_orthrand, pos_orthrand] = computeTime(time_orthrand);
    [times_twosided, neg_twosided, pos_twosided] = computeTime(time_twosided);
    [times_HaTT, neg_HaTT, pos_HaTT] = computeTime(time_HaTT);
    [times_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd] = computeTime(time_HaTT_no_svd);

    subplot(1, 3, 2)

    errorbar(D_set, times_TTrounding,       neg_TTrounding,       pos_TTrounding,       'o-', 'Color', '#3570b6','markersize',6,'linewidth',1.5);

    hold on

    errorbar(D_set, times_randorth,  neg_randorth,  pos_randorth,  '^-', 'Color', '#23a6ba','markersize',6,'linewidth',1.5);
    errorbar(D_set, times_orthrand,  neg_orthrand,  pos_orthrand,  '+-', 'Color', '#b482ba','markersize',6,'linewidth',1.5);
    errorbar(D_set, times_twosided, neg_twosided, pos_twosided, '*-', 'Color', '#c48f00','markersize',6,'linewidth',1.5);
    errorbar(D_set, times_HaTT,       neg_HaTT,       pos_HaTT,       'ks-','markersize',6,'linewidth',1.5);
    errorbar(D_set, times_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd, 'x-', 'Color', '#c45c30','markersize',6,'linewidth',1.5);
    
    hold off
    % title("Alpine-(c)")
    xlabel('Dimension', 'FontSize', 18)
    ylabel('Time (s)', 'FontSize', 18)
    % legend('TTrounding', 'HaTT', 'HaTT-no-SVD', 'randorth', 'orthrand')
    set(gca,'FontSize',14,"FontName", "Times New Roman", 'YScale','log')
    axis square;

        
    % Post-process timings
    
    [speedup_TTrounding, neg_TTrounding, pos_TTrounding] = computeSpeedup(time_TTrounding, time_TTrounding);
    [speedup_randorth, neg_randorth, pos_randorth] = computeSpeedup(time_randorth, time_TTrounding);
    [speedup_orthrand, neg_orthrand, pos_orthrand] = computeSpeedup(time_orthrand, time_TTrounding);
    [speedup_twosided, neg_twosided, pos_twosided] = computeSpeedup(time_twosided, time_TTrounding);
    [speedup_HaTT, neg_HaTT, pos_HaTT] = computeSpeedup(time_HaTT, time_TTrounding);
    [speedup_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd] = computeSpeedup(time_HaTT_no_svd, time_TTrounding);
    
    subplot(1, 3, 3)
    
    errorbar(D_set, speedup_TTrounding,       zeros(length(D_set), 1),       zeros(length(D_set), 1),       'o-', 'Color', '#3570b6','markersize',6,'linewidth',1.5);
    
    hold on
    
    errorbar(D_set, speedup_randorth,  neg_randorth,  pos_randorth,  '^-', 'Color', '#23a6ba','markersize',6,'linewidth',1.5);
    errorbar(D_set, speedup_orthrand,  neg_orthrand,  pos_orthrand,  '+-', 'Color', '#b482ba','markersize',6,'linewidth',1.5);
    errorbar(D_set, speedup_twosided, neg_twosided, pos_twosided, '*-', 'Color', '#c48f00','markersize',6,'linewidth',1.5);
    errorbar(D_set, speedup_HaTT,       neg_HaTT,       pos_HaTT,       'ks-','markersize',6,'linewidth',1.5);
    errorbar(D_set, speedup_HaTT_no_svd, neg_HaTT_no_svd, pos_HaTT_no_svd, 'x-', 'Color', '#c45c30','markersize',6,'linewidth',1.5);
    
    hold off
    % title("Alpine-(b)")
    xlabel('Dimension', 'FontSize', 18)
    ylabel('Speedup', 'FontSize', 18)
    % legend('TTrounding', 'HaTT', 'HaTT-no-SVD', 'randorth', 'orthrand')
    set(gca,'FontSize',14,"FontName", "Times New Roman")
    axis square;

    % % Post-process partical contraction timings
    % [times_HBilinearForm, neg_HBilinearForm, pos_HBilinearForm] = computeTime(time_HBilinearForm);
    % [times_HBilinearForm_no_svd, neg_HBilinearForm_no_svd, pos_HBilinearForm_no_svd] = computeTime(time_HBilinearForm_no_svd);
    % [times_dot, neg_dot, pos_dot] = computeTime(time_dot);
    % 
    % subplot(2, 2, 4)
    % 
    % hold on
    % 
    % errorbar(D_set, times_dot,  neg_dot,  pos_dot,  '^-', 'Color', '#23a6ba','markersize',6,'linewidth',1.5);
    % errorbar(D_set, times_HBilinearForm,       neg_HBilinearForm,       pos_HBilinearForm,       'ks-','markersize',6,'linewidth',1.5);
    % errorbar(D_set, times_HBilinearForm_no_svd, neg_HBilinearForm_no_svd, pos_HBilinearForm_no_svd, 'x-', 'Color', '#c45c30','markersize',6,'linewidth',1.5);
    % 
    % hold off
    % title("(d)")
    % xlabel('Dimension', 'FontSize', 18)
    % ylabel('Time (s)', 'FontSize', 18)
    % legend('dot', 'HBilinearForm-1', 'HBilinearForm-2')
    % set(gca,'FontSize',14,"FontName", "Times New Roman")
    % box on
    % axis square;
end