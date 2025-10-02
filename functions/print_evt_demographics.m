
function print_evt_demographics(PMA_evt_start, PNA_evt_start)
    n = size(PMA_evt_start,2);
    assert(n>0);
    assert(size(PNA_evt_start,2) == n);

    PMA_evt_start_str = sprintf('PMA_evt_start (weeks): \t%.2f (%.2f-%.2f)\t(Missing: %d/%d)', median(PMA_evt_start,'omitnan'), quantile(PMA_evt_start,0.25),quantile(PMA_evt_start,0.75),sum(isnan(PMA_evt_start)),n);
    PNA_evt_start_str = sprintf('PNA_evt_start (weeks): \t%.2f (%.2f-%.2f)\t(Missing: %d/%d)', median(PNA_evt_start,'omitnan'), quantile(PNA_evt_start,0.25),quantile(PNA_evt_start,0.75),sum(isnan(PNA_evt_start)),n);
    fprintf('%s\n%s\n%s\n', PMA_evt_start_str, PNA_evt_start_str);
end