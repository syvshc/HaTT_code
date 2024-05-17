function f = PlotResult_sincos(varargin)
    if nargin == 1
        file = varargin{1};
        name = 'sincos';
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
    
    [err_TTrounding, neg_TTrounding, pos_TTrounding] = computeError(errors_TTrounding);
    [err_randorth, neg_randorth, pos_randorth] = computeError(errors_randorth);
    [err_orthrand, neg_orthrand, pos_orthrand] = computeError(errors_orthrand);
    [err_HaTT1, neg_HaTT1, pos_HaTT1] = computeError(errors_HaTT1);
    [err_HaTT2, neg_HaTT2, pos_HaTT2] = computeError(errors_HaTT2);
    [err_twosided, neg_twosided, pos_twosided] = computeError(errors_twosided);

    subplot(1, 3, 1)


    errorbar(test_ranks, err_TTrounding, neg_TTrounding, pos_TTrounding,     'o', 'Color', '#3570b6','markersize',6,'linewidth',1.5)
    ax = gca;
    ax.YScale = 'log';
    
    hold on
    
    errorbar(test_ranks, err_randorth, neg_randorth, pos_randorth, '^', 'Color', '#23a6ba','markersize',6,'linewidth',1.5)
    errorbar(test_ranks, err_orthrand, neg_orthrand, pos_orthrand, '+', 'Color', '#b482ba','markersize',6,'linewidth',1.5)
    errorbar(test_ranks, err_twosided, neg_twosided, pos_twosided, '*', 'Color', '#c48f00','markersize',6,'linewidth',1.5)
    errorbar(test_ranks, err_HaTT1, neg_HaTT1, pos_HaTT1,     'ks','markersize',6,'linewidth',1.5)
    errorbar(test_ranks, err_HaTT2, neg_HaTT2, pos_HaTT2, 'x', 'Color', '#c45c30','markersize',6,'linewidth',1.5)
    
    hold off
    
    % title('(a)')
    xlabel('Target ranks', 'FontSize', 18)
    ylabel('Relative Error', 'FontSize', 18)
    legend('TT-Rounding', 'RandOrth', 'OrthRand', 'TwoSided', 'HaTT-1', 'HaTT-2')
    set(gca,'FontSize',14,"FontName", "Times New Roman")
    axis square;

    % Post-process timings

    [times_TTrounding, neg_TTrounding, pos_TTrounding] = computeTime(time_TTrounding);
    [times_randorth, neg_randorth, pos_randorth] = computeTime(time_randorth);
    [times_orthrand, neg_orthrand, pos_orthrand] = computeTime(time_orthrand);
    [times_twosided, neg_twosided, pos_twosided] = computeTime(time_twosided);
    [times_HaTT1, neg_HaTT1, pos_HaTT1] = computeTime(time_HaTT1);
    [times_HaTT2, neg_HaTT2, pos_HaTT2] = computeTime(time_HaTT2);

    subplot(1, 3, 2)

    errorbar(test_ranks, times_TTrounding,       neg_TTrounding,       pos_TTrounding,       'o-', 'Color', '#3570b6','markersize',6,'linewidth',1.5);

    hold on

    errorbar(test_ranks, times_randorth,  neg_randorth,  pos_randorth,  '^-', 'Color', '#23a6ba','markersize',6,'linewidth',1.5);
    errorbar(test_ranks, times_orthrand,  neg_orthrand,  pos_orthrand,  '+-', 'Color', '#b482ba','markersize',6,'linewidth',1.5);
    errorbar(test_ranks, times_twosided, neg_twosided, pos_twosided, '*-', 'Color', '#c48f00','markersize',6,'linewidth',1.5);
    errorbar(test_ranks, times_HaTT1,       neg_HaTT1,       pos_HaTT1,       'ks-','markersize',6,'linewidth',1.5);
    errorbar(test_ranks, times_HaTT2, neg_HaTT2, pos_HaTT2, 'x-', 'Color', '#c45c30','markersize',6,'linewidth',1.5);
    
    hold off
    % title("(c)")
    xlabel('Target ranks', 'FontSize', 18)
    ylabel('Time (s)', 'FontSize', 18)
    % legend('TTrounding', 'HaTT', 'HaTT-no-SVD', 'randorth', 'orthrand')
    set(gca,'FontSize',14,"FontName", "Times New Roman")
    axis square;

        % Post-process timings
    
        [speedup_HaTT1, neg_HaTT1, pos_HaTT1] = computeSpeedup(time_HaTT1, time_TTrounding);
        [speedup_HaTT2, neg_HaTT2, pos_HaTT2] = computeSpeedup(time_HaTT2, time_TTrounding);
        [speedup_TTrounding, neg_TTrounding, pos_TTrounding] = computeSpeedup(time_TTrounding, time_TTrounding);
        [speedup_randorth, neg_randorth, pos_randorth] = computeSpeedup(time_randorth, time_TTrounding);
        [speedup_orthrand, neg_orthrand, pos_orthrand] = computeSpeedup(time_orthrand, time_TTrounding);
        [speedup_twosided, neg_twosided, pos_twosided] = computeSpeedup(time_twosided, time_TTrounding);
        
        subplot(1, 3, 3)
        
        errorbar(test_ranks, speedup_TTrounding,       zeros(length(test_ranks), 1),       zeros(length(test_ranks), 1),       'o-', 'Color', '#3570b6','markersize',6,'linewidth',1.5);
        
        hold on
        
        errorbar(test_ranks, speedup_randorth,  neg_randorth,  pos_randorth,  '^-', 'Color', '#23a6ba','markersize',6,'linewidth',1.5);
        errorbar(test_ranks, speedup_orthrand,  neg_orthrand,  pos_orthrand,  '+-', 'Color', '#b482ba','markersize',6,'linewidth',1.5);
        errorbar(test_ranks, speedup_twosided, neg_twosided, pos_twosided, '*-', 'Color', '#c48f00','markersize',6,'linewidth',1.5);
        errorbar(test_ranks, speedup_HaTT1,       neg_HaTT1,       pos_HaTT1,       'ks-','markersize',6,'linewidth',1.5);
        errorbar(test_ranks, speedup_HaTT2, neg_HaTT2, pos_HaTT2, 'x-', 'Color', '#c45c30','markersize',6,'linewidth',1.5);
        
        hold off
        % title("(b)")
        xlabel('Target ranks', 'FontSize', 18)
        ylabel('Speedup', 'FontSize', 18)
        % legend('TTrounding', 'HaTT', 'HaTT-no-SVD', 'randorth', 'orthrand')
        set(gca,'FontSize',14,"FontName", "Times New Roman")
        axis square;
end