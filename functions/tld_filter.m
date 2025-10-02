function tld = tld_filter(tld,a,b)
    % Remove events with event start outside the range [a,b[
    keep_idx = (tld.evt_start_pna_days >=a) & (tld.evt_start_pna_days<b);
    if isfield(tld, 'good_rate_volume')
        keep_idx = keep_idx & tld.good_rate_volume;
    end

    fields_dim1 = {'PMA','BIRTH_weight','hb_pre','hb_post','fiO2_pre','fiO2_post','transfusion_volume',...
        'transfusion_rate','evt_start_pma_days','evt_start_pna_days','TESTOD_MostRecentWeight','hr_mean',...
        'hr_std','sats_mean','sats_std','rr_mean','rr_std','good_rate_volume','sig_quality','pre_post_hb'};
    fields_dim0 = {'t_stop','studyid', 'subjectid'};
    cell_fields = {'BIRTH_Gender','TESTOD_CurrentVentilation'};
    
    for ifieldName = 1:length(fields_dim1)
        fieldName = fields_dim1{ifieldName};
        if isfield(tld,fieldName)
            tld.(fieldName) = tld.(fieldName)(:,keep_idx);
        end
    end
    
    for ifieldName = 1:length(fields_dim0)
        fieldName = fields_dim0{ifieldName};
        if isfield(tld,fieldName)
            tld.(fieldName) = tld.(fieldName)(keep_idx,:);
        end
    end

    for ifieldName = 1:length(cell_fields)
        fieldName = cell_fields{ifieldName};
        if isfield(tld,fieldName)
            tld.(fieldName) = tld.(fieldName)(:,keep_idx);
        end
    end
    tld.counter = sum(keep_idx)+1;
end
