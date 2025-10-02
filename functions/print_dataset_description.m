function print_dataset_description(tld, signals_, dataset_label, plot_dir)
    % Get the index of only one patient ID
    patient_ids_ = tld.subjectid;

    signals = signals_(:,tld.good_rate_volume,:);
    patient_ids = patient_ids_(tld.good_rate_volume);
    PNA_evt_start = tld.evt_start_pna_days(tld.good_rate_volume)/7;
    PMA_evt_start = tld.evt_start_pma_days(tld.good_rate_volume)/7;


    [~,idx] = unique(patient_ids);
    
    
    hb_signals = tld.pre_post_hb(:,tld.good_rate_volume);

    BIRTH_Gender = tld.BIRTH_Gender(idx);
    PMA = tld.PMA(idx) /7;


    BIRTH_weight = tld.BIRTH_weight(idx);
    respirators = tld.TESTOD_CurrentVentilation(tld.good_rate_volume);
    
    a = unique(patient_ids, 'stable');
    count_transfusions = cell2mat(cellfun(@(x) sum(ismember(patient_ids, x)), a,'un',0));
    
    fprintf('Total, [n-evts]=%d, [n-patients]=%d\n', size(signals,2), size(idx,1));
    print_demographics(BIRTH_weight, BIRTH_Gender, PMA);
    print_evt_demographics(PMA_evt_start, PNA_evt_start)
    print_respirator(respirators)
    
    incr_hb = hb_signals(2,:)-hb_signals(1,:);
    incr_hb_str=sprintf('Increment Hb: %.1f (%.1f) (Missing: %d/%d)', mean(hb_signals(2,:) - hb_signals(1,:),'omitnan'), std(hb_signals(2,:)-hb_signals(1,:),'omitnan'),sum(isnan(incr_hb)),size(hb_signals,2) );
    %pre_hb_str=sprintf('Pre Hb: %.1f (%.1f) (Missing: %d/%d)', mean(hb_signals(1,:),'omitnan'), std(hb_signals(1,:),'omitnan'),sum(isnan(hb_signals(1,:))),size(hb_signals,2));
    pre_hb_str = data_description(hb_signals(1,:), "Pre-Hb (g/L)");
    post_hb_str=sprintf('Post Hb: %.1f (%.1f)  (Missing: %d/%d)', mean(hb_signals(2,:),'omitnan'), std(hb_signals(2,:),'omitnan'),sum(isnan(hb_signals(2,:))),size(hb_signals,2));
    
    vol = tld.transfusion_volume(tld.good_rate_volume);
    rates = tld.transfusion_rate(tld.good_rate_volume);
    weights = tld.TESTOD_MostRecentWeight(tld.good_rate_volume);
    n = size(vol,2);
    
    volumes_str = data_description(vol, 'Volumes (ml)'); 
    rates_str = data_description(rates, 'Rates');
    
    volumes_per_kilo_str = data_description(vol./weights, 'Dose (ml/kg)'); 
    rates_per_kg_str = data_description(rates./weights, 'Rates (/kg)');
  
    fig=figure; plot(vol./weights, rates./weights, '.', 'MarkerSize',10)
    xlabel('Dose (ml/kg)'); ylabel('Rates (ml/h/kg)')
    exportgraphics(fig, sprintf('%s/%s_rates_vs_dose.jpg',plot_dir, dataset_label), 'BackgroundColor', 'none','Resolution',300);
    
    transfusion_counts_str = data_description(count_transfusions', '[n-evts] per patients');
    
    fprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n', pre_hb_str, post_hb_str, incr_hb_str,volumes_str,volumes_per_kilo_str,rates_str,rates_per_kg_str,transfusion_counts_str);
    
end