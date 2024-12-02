function f = PlotResult_RDVDAC3(varargin)
    if nargin == 1
        file = varargin{1};
        name = 'RDVDAC3';
    else
        file = varargin{1};
        name = varargin{2};
    end
    % load(file, "errors*", "d_set", "errors*", "time*", "energy*");
    load(file)
    %% Plot results
    
    f = figure('Name', name);
    % f.Position(1:2) = [0,1050];
    % f.Position(3:4) = [1050, 700];
    f.Position = [1,49,1440,781.5];
    % Post-process errors

    errors_TTrounding = errors_TTrounding * ones(1, 7);
    time_pcg = time_pcg * ones(1, 7);
    
    [err_TTrounding, neg_TTrounding, pos_TTrounding] = computeError(errors_TTrounding);
    [err_randorth, neg_randorth, pos_randorth] = computeError(errors_randorth);
    [err_orthrand, neg_orthrand, pos_orthrand] = computeError(errors_orthrand);
    [err_twosided, neg_twosided, pos_twosided] = computeError(errors_twosided);
    [err_HaTT1, neg_HaTT1, pos_HaTT1] = computeError(errors_HaTT1);
    [err_HaTT2, neg_HaTT2, pos_HaTT2] = computeError(errors_HaTT2);

    subplot(1, 3, 1)
    
    errorbar(d_set, err_TTrounding, neg_TTrounding, pos_TTrounding, 'o', 'Color', '#3570b6','markersize',6,'linewidth',1.5)
    % ax = gca;
    % ax.YScale = 'log';
    
    hold on
    
    errorbar(d_set, err_randorth, neg_randorth, pos_randorth, '^', 'Color', '#23a6ba','markersize',6,'linewidth',1.5)
    errorbar(d_set, err_orthrand, neg_orthrand, pos_orthrand, '+', 'Color', '#b482ba','markersize',6,'linewidth',1.5)
    errorbar(d_set, err_twosided, neg_twosided, pos_twosided, '*', 'Color', '#c48f00','markersize',6,'linewidth',1.5)
    errorbar(d_set, err_HaTT1, neg_HaTT1, pos_HaTT1,     'ks','markersize',6,'linewidth',1.5)
    errorbar(d_set, err_HaTT2, neg_HaTT2, pos_HaTT2, 'x', 'Color', '#c45c30','markersize',6,'linewidth',1.5)
    
    hold off
    
    % title('(a)')
    xlabel('d''s set', 'FontSize', 18)
    ylabel('Relative Error', 'FontSize', 18)
    legend('TT-Rounding', 'RandOrth', 'OrthRand', 'TwoSided', 'HaTT-1', 'HaTT-2')
    set(gca,'FontSize',14,"FontName", "Times New Roman")
    axis square;

    % Post-process timings
    [times_pcg, neg_pcg, pos_pcg] = computeTime(time_pcg);
    [times_TTrounding, neg_TTrounding, pos_TTrounding] = computeTime(time_TTrounding);
    [times_randorth, neg_randorth, pos_randorth] = computeTime(time_randorth);
    [times_orthrand, neg_orthrand, pos_orthrand] = computeTime(time_orthrand);
    [times_twosided, neg_twosided, pos_twosided] = computeTime(time_twosided);
    [times_HaTT1, neg_HaTT1, pos_HaTT1] = computeTime(time_HaTT1);
    [times_HaTT2, neg_HaTT2, pos_HaTT2] = computeTime(time_HaTT2);

    subplot(1, 3, 2)

    errorbar(d_set, times_TTrounding,       neg_TTrounding,       pos_TTrounding,       'o-', 'Color', '#3570b6','markersize',6,'linewidth',1.5);

    hold on
    errorbar(d_set, times_pcg,  neg_pcg,  pos_pcg,  'o-', 'Color', '#23a6ba','markersize',6,'linewidth',1.5);
    errorbar(d_set, times_randorth,  neg_randorth,  pos_randorth,  '^-', 'Color', '#23a6ba','markersize',6,'linewidth',1.5);
    errorbar(d_set, times_orthrand,  neg_orthrand,  pos_orthrand,  '+-', 'Color', '#b482ba','markersize',6,'linewidth',1.5);
    errorbar(d_set, times_twosided, neg_twosided, pos_twosided, '*-', 'Color', '#c48f00','markersize',6,'linewidth',1.5);
    errorbar(d_set, times_HaTT1,       neg_HaTT1,       pos_HaTT1,       'ks-','markersize',6,'linewidth',1.5);
    errorbar(d_set, times_HaTT2, neg_HaTT2, pos_HaTT2, 'x-', 'Color', '#c45c30','markersize',6,'linewidth',1.5);
    
    hold off
    % title("(c)")
    xlabel('d''s set', 'FontSize', 18)
    ylabel('Time (s)', 'FontSize', 18)
    % legend('TTrounding', 'HaTT', 'HaTT-no-SVD', 'randorth', 'orthrand')
    set(gca,'FontSize',14,"FontName", "Times New Roman", 'YScale','log')
    axis square;

        
    % Post-process timings
    
    % Post-process timings
    % [times_pcg, neg_pcg, pos_pcg] = computeTime(time_pcg);
    [round_times_TTrounding, neg_TTrounding, pos_TTrounding] = computeTime(round_time_TTrounding);
    [round_times_randorth, neg_randorth, pos_randorth] = computeTime(round_time_randorth);
    [round_times_orthrand, neg_orthrand, pos_orthrand] = computeTime(round_time_orthrand);
    [round_times_twosided, neg_twosided, pos_twosided] = computeTime(round_time_twosided);
    [round_times_HaTT1, neg_HaTT1, pos_HaTT1] = computeTime(round_time_HaTT1);
    [round_times_HaTT2, neg_HaTT2, pos_HaTT2] = computeTime(round_time_HaTT2);

    subplot(1, 3, 3)

    errorbar(d_set, round_times_TTrounding,       neg_TTrounding,       pos_TTrounding,       'o-', 'Color', '#3570b6','markersize',6,'linewidth',1.5);

    hold on
    % errorbar(d_set, round_times_pcg,  neg_pcg,  pos_pcg,  'o-', 'Color', '#23a6ba','markersize',6,'linewidth',1.5);
    errorbar(d_set, round_times_randorth,  neg_randorth,  pos_randorth,  '^-', 'Color', '#23a6ba','markersize',6,'linewidth',1.5);
    errorbar(d_set, round_times_orthrand,  neg_orthrand,  pos_orthrand,  '+-', 'Color', '#b482ba','markersize',6,'linewidth',1.5);
    errorbar(d_set, round_times_twosided, neg_twosided, pos_twosided, '*-', 'Color', '#c48f00','markersize',6,'linewidth',1.5);
    errorbar(d_set, round_times_HaTT1,       neg_HaTT1,       pos_HaTT1,       'ks-','markersize',6,'linewidth',1.5);
    errorbar(d_set, round_times_HaTT2, neg_HaTT2, pos_HaTT2, 'x-', 'Color', '#c45c30','markersize',6,'linewidth',1.5);
    
    hold off
    % title("(b)")
    xlabel('d''s set', 'FontSize', 18)
    ylabel('Speedup', 'FontSize', 18)
    % legend('TTrounding', 'HaTT', 'HaTT-no-SVD', 'randorth', 'orthrand')
    set(gca,'FontSize',14,"FontName", "Times New Roman")
    axis square;
end