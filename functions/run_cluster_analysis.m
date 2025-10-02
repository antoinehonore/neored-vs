function run_cluster_analysis(tld,t_overlap,t_limits)
    pre = find(tld.t_average<0);
    post = find(tld.t_average>0);
    sig_quality = tld.sig_quality;
    
    if strcmp(tld.dataset_label, 'Stockholm')
        idx_include_hr = sig_quality(1,:);
        idx_include_sats = sig_quality(2,:);
        idx_include_rr = sig_quality(3,:);
    elseif strcmp(tld.dataset_label, 'Imperial')
        idx_include_hr = sig_quality(1,:);
        idx_include_sats = sig_quality(2,:);
        idx_include_rr = sig_quality(3,:);
    else
        idx_include_hr = idx_include;
        idx_include_sats = idx_include;
        idx_include_rr = idx_include;
    end

    % dataPrehrmean=tld.hr_mean(pre,idx_include_hr);
    % dataPosthrmean=tld.hr_mean(post,idx_include_hr);
    % disp('cluster analysis HR mean') 
    % montecarlo_fornewdata_differentsizes(dataPrehrmean', dataPosthrmean', 1/t_overlap, t_limits(1)+1, t_limits(2));
    % 
    % dataPrehrsd=tld.hr_std(pre,idx_include_hr);
    % dataPosthrsd=tld.hr_std(post,idx_include_hr);
    % disp('cluster analysis HR SD') 
    % montecarlo_fornewdata_differentsizes(dataPrehrsd', dataPosthrsd', 1/t_overlap, t_limits(1)+1, t_limits(2));
    % 
    % dataPresatsmean=tld.sats_mean(pre,idx_include_sats);
    % dataPostsatsmean=tld.sats_mean(post,idx_include_sats);
    % disp('cluster analysis Sats mean') 
    % montecarlo_fornewdata_differentsizes(dataPresatsmean', dataPostsatsmean', 1/t_overlap, t_limits(1)+1, t_limits(2));
    % 
    % dataPresatssd=tld.sats_std(pre,idx_include_sats);
    % dataPostsatssd=tld.sats_std(post,idx_include_sats);
    % disp('cluster analysis Sats SD') 
    % montecarlo_fornewdata_differentsizes(dataPresatssd', dataPostsatssd', 1/t_overlap, t_limits(1)+1, t_limits(2));
    % 
    % dataPrerrmean=tld.rr_mean(pre,idx_include_rr);
    % dataPostrrmean=tld.rr_mean(post,idx_include_rr);
    % disp('cluster analysis RR mean')
    % montecarlo_fornewdata_differentsizes(dataPrerrmean', dataPostrrmean', 1/t_overlap, t_limits(1)+1, t_limits(2));
    % 
    % dataPrerrsd=tld.rr_std(pre,idx_include_rr);
    % dataPostrrsd=tld.rr_std(post,idx_include_rr);
    % disp('cluster analysis RR SD')
    % montecarlo_fornewdata_differentsizes(dataPrerrsd', dataPostrrsd', 1/t_overlap, t_limits(1)+1, t_limits(2));

    diary('cluster_analysis_results.txt');

    results = {};

    % --- HR mean ---
    dataPrehrmean = tld.hr_mean(pre,idx_include_hr);
    dataPosthrmean = tld.hr_mean(post,idx_include_hr);
    disp('cluster analysis HR mean')
    resHRmean = montecarlo_fornewdata_differentsizes(dataPrehrmean', dataPosthrmean', 1/t_overlap, t_limits(1)+1, t_limits(2));
    results(end+1,:) = {'HR mean', resHRmean};
    class(resHRmean)

    % --- HR SD ---
    dataPrehrsd = tld.hr_std(pre,idx_include_hr);
    dataPosthrsd = tld.hr_std(post,idx_include_hr);
    disp('cluster analysis HR SD')
    resHRsd = montecarlo_fornewdata_differentsizes(dataPrehrsd', dataPosthrsd', 1/t_overlap, t_limits(1)+1, t_limits(2));
    results(end+1,:) = {'HR SD', resHRsd};

    % --- Sats mean ---
    dataPresatsmean = tld.sats_mean(pre,idx_include_sats);
    dataPostsatsmean = tld.sats_mean(post,idx_include_sats);
    disp('cluster analysis Sats mean')
    resSatsMean = montecarlo_fornewdata_differentsizes(dataPresatsmean', dataPostsatsmean', 1/t_overlap, t_limits(1)+1, t_limits(2));
    results(end+1,:) = {'Sats mean', resSatsMean};

    % --- Sats SD ---
    dataPresatssd = tld.sats_std(pre,idx_include_sats);
    dataPostsatssd = tld.sats_std(post,idx_include_sats);
    disp('cluster analysis Sats SD')
    resSatsSD = montecarlo_fornewdata_differentsizes(dataPresatssd', dataPostsatssd', 1/t_overlap, t_limits(1)+1, t_limits(2));
    results(end+1,:) = {'Sats SD', resSatsSD};

    % --- RR mean ---
    dataPrerrmean = tld.rr_mean(pre,idx_include_rr);
    dataPostrrmean = tld.rr_mean(post,idx_include_rr);
    disp('cluster analysis RR mean')
    resRRmean = montecarlo_fornewdata_differentsizes(dataPrerrmean', dataPostrrmean', 1/t_overlap, t_limits(1)+1, t_limits(2));
    results(end+1,:) = {'RR mean', resRRmean};

    % --- RR SD ---
    dataPrerrsd = tld.rr_std(pre,idx_include_rr);
    dataPostrrsd = tld.rr_std(post,idx_include_rr);
    disp('cluster analysis RR SD')
    resRRsd = montecarlo_fornewdata_differentsizes(dataPrerrsd', dataPostrrsd', 1/t_overlap, t_limits(1)+1, t_limits(2));
    results(end+1,:) = {'RR SD', resRRsd};

    % Stop capturing diary
    diary off

    %% Collect results into a CSV
    measures = {'HR mean','HR SD','Sats mean','Sats SD','RR mean','RR SD'};
    resTables = {resHRmean,resHRsd,resSatsMean,resSatsSD,resRRmean,resRRsd};

    allTables = cell(numel(resTables),1);

    for i = 1:numel(resTables)
        thisTable = resTables{i};
        if istable(thisTable)
            % Add Measure column
            thisTable = addvars(thisTable, repmat(string(measures{i}),height(thisTable),1), ...
                'Before',1, 'NewVariableNames','Measure');
            allTables{i} = thisTable;
        else
            % In case the function did not return a table
            warning('Result for %s is not a table, skipping', measures{i});
        end
    end

    % Combine into one big table
    finalTable = vertcat(allTables{:});

    whos finalTable %for debugging
    class (finalTable) %for debugging
    head(finalTable) %for debugging


    % Save CSV
 writetable(finalTable,'cluster_analysis_results.csv');

end

