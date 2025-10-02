function tld=substract_mean_baseline(tld,t_baseline,var_labels)
    % mean subtract the moving average to mean of t_baseline
    for v = 1 : numel(var_labels)
        idx_baseline = tld.t_average > t_baseline(1) & tld.t_average < t_baseline(2);
        tld.(var_labels{v}) = tld.(var_labels{v}) - mean(tld.(var_labels{v})(idx_baseline, :), 'omitnan');
    end
end

