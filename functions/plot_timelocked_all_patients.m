function plot_timelocked_all_patients(tld,t_limits, var_labels, y_names, plot_dir)

    dataset_label = tld.dataset_label;

    % plot mean response over all transfusions
    fig = figure;
    po = get(gcf, 'position');
    set(gcf, 'position', [po(1:2), 1600, 1000], 'name', 'timelocked_responses');
    sig_quality = tld.sig_quality;
    
    if isfield(tld, 'ibi_locked'); y_names = [y_names, {'Resp rate [breaths/min]', 'Std IBI [sec]', 'Apnoea rate 5 sec [times/hour]', 'Apnoea rate 10 sec [times/hour]'}]; end
    %if isfield(tld, 'fiO2_mean'); y_names = [y_names, {'Mean fiO2 (%)','Mean Hb', 'Std fiO2 (%)','Std Hb'}]; end
    
    for v = 1 : numel(var_labels)
    
        % plot mean and standard deviations
        subplot(2, numel(y_names) / 2, v);
        
        idx_include = find(sig_quality(v, :));
        
        m = mean(tld.(var_labels{v})(:, idx_include), 2, 'omitnan');
        s = std(tld.(var_labels{v})(:, idx_include), [], 2, 'omitnan');
        plot(tld.t_average, m, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
        hold on;
        patch([tld.t_average; tld.t_average(end : -1 : 1)], [m + s; m(end : -1 : 1) - s(end : -1 : 1)], [0.5, 0.5, 0.5], 'FaceAlpha', 0.5, 'EdgeColor', [1, 1, 1]);
        plot([0, 0], [min(m - s), max(m + s)], 'k--', 'LineWidth', 2);
        plot(t_limits, [0, 0], 'k--', 'LineWidth', 1);
        %plot([mean(tld.t_stop(idx_include), 'omitnan'), mean(tld.t_stop(idx_include), 'omitnan')], [min(m - s), max(m + s)], 'k--', 'LineWidth', 2);
        %instead of previous version (above) make sure to plot the mean stop
        %time and deal with missing data by inserting standard dotted instead of dashed line at +4h 
        stop_time = mean(tld.t_stop(idx_include), 'omitnan');
        ymin = min(m - s);
        ymax = max(m + s);
    
        if ~isnan(stop_time) && isfinite(stop_time)
            plot([stop_time, stop_time], [ymin, ymax], 'k--', 'LineWidth', 2);
        else
            % Fallback: vertical dotted line at x = 4
            plot([4, 4], [ymin, ymax], 'k:', 'LineWidth', 1);
        end
        title(sprintf('%s', var_labels{v}), 'Interpreter', 'None');
        xlabel('Time [hours]');
        ylabel(y_names{v});
        xlim([tld.t_average([1, end])])
        ylim([ymin, ymax])
    end
    figure(fig);
    orient(fig, 'landscape');
    exportgraphics(fig, sprintf('%s/timelocked_data_all_patients_%s.jpg',plot_dir,dataset_label), 'BackgroundColor', 'none','Resolution',300);
    print(fig, '-dpdf', '-bestfit', sprintf('%s/timelocked_data_all_patients_%s.pdf', plot_dir, dataset_label));
    savefig(fig, sprintf('%s/timelocked_data_all_patients_%s.fig',plot_dir,dataset_label) );


end

