function [signals, colors,signal_increase,signal_decrease,var_labels,tld]=plot_vital_signs(files, t_limits, t_window, t_overlap, t_baseline, sig_check, dataset_label, std_threshold, plot_dir, ...
    sig_quality_fname,transfusion_rawdata_fname,evt_start_days_lower,evt_start_days_upper,load_other_centers,cluster_analysis)
% - plot_vital_signs(files, t_window, t_overlap, sig_check, dataset_label,
% std_threshold) -- plots the vital signs responses after time-locking to 
% transfusion events.
% 
% 
% Input
% - files           - File paths to be analysed. Paths should be defined as
%                     directory names in cell structure.
% - t_limits        - Time limits in hours defined as double variable
% - t_window        - Window size in hours defined as double variable.
% - t_overlap       - Overlap between consecutive windows in hours defined 
%                     as double.
% - t_baseline      - Time limits in hours used for baseline correction
%                     (i.e., mean subtraction) defined as double variable.
% - sig_check       - Boolean specifying if vital signs time series should
%                     be visualise and checked.
% - dataset_label   - 'Berlin', 'Oxford', or 'Stockholm'. Will be used to
%                     reformat data structure if necessary.
% - std_threhsold   - Threshold used to define (non-)responders;
%                     it should be exceeded post-transfusion start to be
%                     classified as responder.
%
%
% Example usage
% - define_subjects_Oxford (example script for data structure).
% - t_window = 1; % hours
% - t_overlap = 0.5; % hours
% - t_limits = [-6, 12]; % hours
% - t_baseline = [-6, 0]; % hours
% - check_sig = true;
% - dataset_label = 'Oxford';
% - std_threshold = 1.5;
% - plot_vital_signs(files, t_limits, t_window, t_overlap, t_baseline, ...
%       check_sig, dataset_label, std_threshold)
%
%
% See also MOVING_AVERAGE, TIME_LOCK_DATA, CHECK_SIGNALS, FIND_RESPONDERS
% ________________________________________________________________________
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html
%
%                                           (c) Coen Zandvoort, Caroline Hartley 2024
% ________________________________________________________________________

% check input and set defaults if necessary
if nargin < 2, t_limits = [-6, 12]; end
if nargin < 3, t_window = 1; end
if nargin < 4, t_overlap = 0.5; end
if nargin < 5, sig_check = true; end
if nargin < 6, dataset_label = 'Oxford'; end
if nargin < 7, std_threshold = 1.5; end


plot_subject_specific=false;
auto_sig_check=1;

% check the availability of the functions folder
if exist(transfusion_rawdata_fname)
    load(transfusion_rawdata_fname)
else
    % pre-allocate some struct and fields
    tld = struct;
    tld.counter = 1;
    tld_nocorrection = struct;
    tld_nocorrection.counter=1;
    tld.BIRTH_Gender = {};
    tld.PMA = [];
    tld.BIRTH_weight = [];
    tld.hb_pre = [];
    tld.hb_post = [];
    tld.fiO2_pre = [];
    tld.fiO2_post = [];
    tld.transfusion_volume = [];
    tld.transfusion_rate = [];

    tld.evt_start_pma_days = [];
    tld.evt_start_pna_days = [];

    tld.TESTOD_CurrentVentilation={};
    tld.TESTOD_MostRecentWeight=[];
    
    toremove=[]; subjectid=[]; studyid=[];
    
    var_labels = {'hr_mean', 'sats_mean', 'rr_mean', 'hr_std', 'sats_std', 'rr_std'};
    
    if strcmp(dataset_label, 'Stockholm')
        subjectid= cell(size(files, 1),1);
    end
    
    for f = 1 : size(files, 1)
        
        % load vital signs and inter-breath interval data
        %fprintf('reading:\t"%s"\n', files{f, 1})
        vs = load(files{f, 1});

        if size(files, 2) > 1
            if ~isempty(files{f, 2})
                % load inter-breath intervals
                vs.ibi = load(files{f, 2});
            end
        end
        
        % Berlin and Stockholm data structure should be reformatted to match 
        % the Oxford data
        vs.dataset_label = dataset_label;
        
        if strcmp(dataset_label, 'Stockholm')
            [vs,studyid,subjectid{f}] = standardize_stockholm(vs,studyid,files,f);
        end
        
        if strcmp(dataset_label, 'Berlin')
            % Here, we may have to add some extra lines of code to make the
            % Berlin and Stockholm data similar to the Oxford data. 
            vs.time_events_vital_signs = vs.time_events_vital_signs;
            vs.events_vital_signs = {'rbctransfusion_start','rbctransfusion_stop'};
            vs.HR = vs.HR;
            vs.sats = vs.sats;
            vs.RR = vs.RR;
            vs.BP = vs.BP;
            vs.fiO2 = vs.fiO2;
        end
    
        if strcmp(dataset_label, 'Imperial')
            % Here, we may have to add some extra lines of code to make the
            % Berlin and Stockholm data similar to the Oxford data. 
            vs.time_events_vital_signs = vs.time_events_vital_signs;
            vs.events_vital_signs = {'rbctransfusion_start','rbctransfusion_stop'};
            vs.HR = vs.HR;
            vs.sats = vs.sats_max; %change sats measure here, toggle between median, max and just 'sats' which is equal to median
            % vs.sats_max = vs.sats_max;
            % vs.sats_median = vs.sats_median;
            vs.rate=vs.volume/vs.duration;
            vs.birthweight = vs.BIRTH_weight;
            vs.PMA = vs.PMA;
            vs.pna = vs.pna;
            vs.sex = vs.BIRTH_Gender;
            vs.PMA_at_start = vs.PMA_at_start;
            vs.evt_sum_dose = vs.volume;
            vs.hb_post = vs.post_Hb;
            vs.hb_pre = vs.pre_Hb;
            vs.RR = vs.RR;
            vs.BP = vs.BP;
            vs.fiO2 = vs.fiO2;
        end
        
        % find transfusion starts
        % find transfusion starts - will skip data set if it can't find
        % transfusion start
        if ~isfield(vs, 'events_vital_signs'); disp('no data'); toremove=[toremove;f]; continue; end
        
        if ~any(find(contains(vs.events_vital_signs, 'transfusion') & ...
            contains(vs.events_vital_signs, 'start'))); disp('no transfusion marker'); toremove=[toremove;f]; continue; end
        
        if strcmp(dataset_label, 'Berlin')
            fs=0.0011;
        end
    
        if strcmp(dataset_label, 'Oxford')
            fs = 0.9766;
        end
        
        if  strcmp(dataset_label, 'Stockholm')
            fs=1;
        end

        if strcmp(dataset_label, 'Imperial')
            fs=0.016;
        end
        
        % time-lock data
        tld = time_lock_data(vs, tld, t_limits, t_baseline, fs);

        if any(strcmp(dataset_label, {'Oxford', 'Stockholm','Berlin'})) %for high-frequency centres also without correction to use identify_desatbradytachy
            tld_nocorrection = time_lock_data_nocorrection(vs, tld_nocorrection, t_limits, fs);
        end
        
        if strcmp(dataset_label, 'Oxford')
            n = length(find(contains(vs.events_vital_signs, 'transfusion') & contains(vs.events_vital_signs, 'start')));
            for i = 1 : n
                subjectid = [subjectid;files{f}(31:37)];
                studyid = [studyid;files{f}(39:47)];
            end
        end
        if strcmp(dataset_label, 'Imperial')
            n = length(find(contains(vs.events_vital_signs, 'transfusion') & contains(vs.events_vital_signs, 'start')));
            for i = 1 : n
                idx_slash = strfind(files{f}, '/');
                file_name = files{f}(idx_slash(end) + 1 : end);
                idx_underscore = strfind(file_name, 'Imperial_');
                idx_underscore_end = strfind(file_name, '_R_data_');
                file_name = file_name(10 : idx_underscore_end(1) -1);
    
                subjectid = [subjectid; {file_name(1 : 4)}];
                studyid = [studyid; {file_name}];
            end

        end
    end
   
    % get moving average
    tld = moving_average(tld, t_window, t_overlap);
    
    % check vital signs signals visually (note that if sig_check is set to
    % false all data will be included). 
    if isfield(tld, 'ibi_locked'); var_labels = [var_labels, {'resp_rate', 'ibi_std', 'apnoea_rate_5_sec', 'apnoea_rate_10_sec'}]; end
    %if isfield(tld, 'fiO2_mean'); var_labels = [var_labels, {'fiO2_mean','Hb_mean', 'fiO2_std','Hb_std'}]; end
    
    sig_quality = true(numel(var_labels), tld.counter - 1);

    if sig_check
        % extract signal quality (for now, the code will check if the
        % sig_quality.mat file is available. If this is the case, the user
        % doesn't get the option to re do the signal selection. Delete the file
        % sig_quality.mat if you would like to re do this selection).
        if ~exist(sig_quality_fname, 'file')

            sig_quality = true(numel(var_labels), tld.counter - 1);  % Preallocate as requested

            if auto_sig_check
                %automatic checking
                %must have at least 3 hours of baseline in the 6 hours before start and 6 hours post - can
                %change this in next 2 lines
                pre_min_ind=3*t_window/t_overlap;
                post_min_ind=6*t_window/t_overlap;
                pre_ind=intersect(find(tld.t_average<0),find(tld.t_average>=-6));
                post_ind=find(tld.t_average>=0);
                for v = 1 : numel(var_labels)
                    x = tld.(sprintf('%s', var_labels{v}));
                    for u=1:size(x,2)
                        if length(find(~isnan(x(pre_ind,u))))>=pre_min_ind && length(find(~isnan(x(post_ind,u))))>=post_min_ind
                            sig_quality(v, u)=true;
                        else
                            sig_quality(v,u)=false;
                        end
                    end
                end
            else
                % Define variable pairs: mean and std for each type
                paired_vars = {'hr_mean', 'hr_std'; 'sats_mean', 'sats_std'; 'rr_mean', 'rr_std'};

                for p = 1 : size(paired_vars, 1)
                    x1 = tld.(paired_vars{p, 1});
                    x2 = tld.(paired_vars{p, 2});
                    x = cat(3, x1, x2);  % Combine into 3D array
                    time = tld.t_average;
                    events = tld.t_stop;

                    fprintf('\nassessing "%s" and "%s"\n', paired_vars{p, 1}, paired_vars{p, 2});

                    if ~strcmp(dataset_label, 'Stockholm')
                        [found1, found2] = check_signals_dual(x, time, events, t_limits);
                    else
                        found1 = sum(isnan(x1), 1) == 0;
                        found2 = sum(isnan(x2), 1) == 0;
                    end

                    % Store results in preallocated sig_quality
                    idx1 = find(strcmp(var_labels, paired_vars{p, 1}));
                    idx2 = find(strcmp(var_labels, paired_vars{p, 2}));
                    sig_quality(idx1, :) = found1;
                    sig_quality(idx2, :) = found2;
                end
            end

            % Save output
            save(sig_quality_fname, 'sig_quality');
        end

end
    close
    load(sig_quality_fname, 'sig_quality')
    
    alpha = 0.2;
    timeline = tld.t_average;
    %clear std_threshold
    save(transfusion_rawdata_fname)
end

%print(sprintf('../results/timelocked_data_all_patients_%s_2.jpg',dataset_label), '-djpg', '-bestfit')

%all_markers = {'s','*','o','p'};
rates_per_kilo = tld.transfusion_rate./tld.TESTOD_MostRecentWeight;
vol_per_kilo = tld.transfusion_volume./tld.TESTOD_MostRecentWeight;

tld.good_rate_volume = ((rates_per_kilo >= 3) & vol_per_kilo >= 8)|(isnan(rates_per_kilo)|isnan(vol_per_kilo));

sig_quality(:, ~tld.good_rate_volume) = 0;


tld.sig_quality = sig_quality;
tld.pre_post_hb = nan(2,tld.counter-1);
tld.pre_post_hb(1,:)=tld.hb_pre;
tld.pre_post_hb(2,:)=tld.hb_post;
tld.dataset_label = dataset_label;
tld.subjectid = subjectid;
tld.studyid = studyid;

% identify bradycardia, tachycardia, desats
if any(strcmp(dataset_label, {'Oxford', 'Stockholm','Berlin'})) %for high-frequency centres
    time_events = identify_desatbradytachy(tld_nocorrection);
    % Hr, Sats, rr
    hr_quality = tld.sig_quality(1, :);
    sats_quality = tld.sig_quality(2, :);

    time_events.desat = time_events.desat(sats_quality);
    time_events.tachy = time_events.tachy(hr_quality);
    time_events.brady = time_events.brady(hr_quality);
    save(sprintf('%s/time_events_%s.mat',plot_dir,dataset_label), 'time_events');
end

y_names = {'Mean HR [bpm]', 'Mean Sats [%]', 'Mean RR [breaths/min]', 'Std HR [bpm]', 'Std Sats [%]', 'Std RR [breaths/min]'};


%Cluster analysis

if cluster_analysis
   run_cluster_analysis(tld, t_overlap, t_limits);
end
tld = rmfield(tld, 'hr');
tld = rmfield(tld, 'sats');
tld = rmfield(tld, 'rr');

save(sprintf('%s/tld_%s.mat',plot_dir,dataset_label), "tld")

% find responders
if 1
    %figure;
    %plot(mean(tld.hr_std,2,'omitnan')); 
    
    tld = substract_mean_baseline(tld,t_baseline,var_labels);
    
    %figure;
    %plot(mean(tld.hr_std,2,'omitnan'));
    
    plot_timelocked_all_patients(tld, t_limits, var_labels, y_names, plot_dir);
    close all;
    tld = tld_filter(tld,evt_start_days_lower, evt_start_days_upper);

    if load_other_centers
        tld_oxford = load('../others/tld_Oxford.mat'); 
        tld_oxford = tld_oxford.tld;

        tld_oxford.dataset_label="Oxford";
        % I don't have the exact PNA for oxford: all the babies are >14
        % days, assume 15
        tld_oxford.evt_start_pna_days = 15*ones(1,tld_oxford.counter-1);

        tld_oxford = substract_mean_baseline(tld_oxford, t_baseline, var_labels);
        plot_timelocked_all_patients(tld_oxford, t_limits, var_labels, y_names, plot_dir);
        
        tld_oxford = tld_filter(tld_oxford, evt_start_days_lower, evt_start_days_upper);

        tld_Imperial = load('../others/tld_Imperial.mat'); 
        tld_Imperial = tld_Imperial.tld;
        tld_Imperial.dataset_label="Imperial";
        tld_Imperial = substract_mean_baseline(tld_Imperial,t_baseline,var_labels);
        plot_timelocked_all_patients(tld_Imperial, t_limits, var_labels, y_names, plot_dir);
        tld_Imperial = tld_filter(tld_Imperial,evt_start_days_lower, evt_start_days_upper);

        %tld = tld_Imperial;
        tld = concatenate_tld(tld, tld_oxford);
        tld = concatenate_tld(tld, tld_Imperial);
    end

    [signals, colors, all_colors, signal_increase, signal_decrease,  all_colors_names] = find_responders(tld, std_threshold, 6, plot_subject_specific); 
    
    %find subject labels for Oxford data
    if strcmp(dataset_label, 'Oxford')
        all_unique_patid = unique(tld.subjectid);
        patient_ids_integer = str2num(tld.subjectid(:,end-1:end));

    end

    if strcmp(dataset_label, 'Imperial')
        all_unique_patid=unique(tld.subjectid);
        patient_ids_integer = 1 : numel(tld.subjectid);
    end
    
    if strcmp(dataset_label, 'Stockholm')
        % Convert the list of hex IDs to integer ids
        
        all_unique_patid = unique(tld.subjectid);
        patient_ids_integer = zeros(size(tld.subjectid,1),1);
        for ipat = 1:size(all_unique_patid,1)
            patient_ids_integer(find(contains(tld.subjectid, all_unique_patid{ipat}))) = ipat;
        end
        
        %print_dataset_description(tld, signals, dataset_label, plot_dir)
    end
    
    % Average lines per response patient groups
    for isig = 1:numel(var_labels)
        fig = figure;
        
        % Good quality recordings
        idx_good_quality = tld.sig_quality(isig, :) == 1;
        if strcmp(dataset_label, 'Oxford') || strcmp(dataset_label, 'Stockholm') || strcmp(dataset_label, 'Imperial')
            % Count unique patients with good quality recordings of the VS
            npatsig = size(unique(patient_ids_integer(idx_good_quality)), 1);
            idx_full_hb = filter_full_hb(tld.pre_post_hb);
            
            full_hb = tld.pre_post_hb(:, idx_good_quality & idx_full_hb);
            %full_hb = hb(:, sum(isnan(hb),1)==0);
            
            patid_included = patient_ids_integer(idx_good_quality);
            [~, idx] = unique(patid_included);
            unique_patid_sig_hb = patid_included(idx);
            fprintf('===============================\n\n');
            fprintf('Variable=%s, n-patients=%d, n-evt=%d, n-fullHB=%d\n', ...
                    var_labels{isig}, npatsig, sum(idx_good_quality), size(full_hb,2));
            
            %print_demographics(tld.BIRTH_weight(unique_patid_sig_hb), tld.BIRTH_Gender(unique_patid_sig_hb), tld.PMA(unique_patid_sig_hb))
        
        end
        % for each types of responses, excluding the "green" for readability (both
        % increase and decrease)
        for icolor = 1:size(all_colors,2)
            
            idx_sigcolor = (idx_good_quality) & (colors(isig,:)==icolor);
            %hb = tld.pre_post_hb(:, idx_sigcolor);
            %full_hb = hb(:, filter_full_hb(hb)); %hb(:, sum(isnan(hb),1)==0);
            
            displayname = '# Transfusions=%d\n# Patients=%d\n';

            if strcmp(dataset_label, 'Oxford') || strcmp(dataset_label, 'Stockholm')
                patid_included = patient_ids_integer(idx_sigcolor);
                [~, idx] = unique(patid_included);
                unique_patid_sigcolor = patid_included(idx);
                
                npatsigcolor = size(idx, 1);
                
                displayname = sprintf(displayname, sum(idx_sigcolor), npatsigcolor);
            end
            % fprintf('\tcolor=%s\n\t%s\n', all_colors_names{icolor}, displayname);
            
            if strcmp(dataset_label, 'Imperial')
                patid_included = patient_ids_integer(idx_sigcolor);    
               % patid_included = subjectid(idx_sigcolor);
                [~, idx] = unique(patid_included);
                unique_patid_sigcolor = patid_included(idx);
                    
                npatsigcolor = size(idx,1);
    
                displayname = sprintf(displayname, sum(idx_sigcolor), npatsigcolor );

            end
            %fprintf('color=%s\n%s\n', all_colors_names{icolor},
            %displayname); @Caroline
            %fprintf('color=%s\n%s\n', all_colors{icolor}, displayname);
            
            % set(findall(0, 'type', 'axes'), 'FontName', 'Times', 'Fontsize', 16, 'TickDir', 'out', 'box', 'off', 'linewidth', 2, 'ticklength', [0.01, 0.01])
                
           %print_demographics(tld.BIRTH_weight(unique_patid_sigcolor),tld.BIRTH_Gender(unique_patid_sigcolor),tld.PMA(unique_patid_sigcolor))
            
           if ~strcmp(all_colors{icolor},'green')
                data_sigcolor = squeeze(signals(isig, idx_sigcolor,:));
                if size(data_sigcolor, 2)==1 %this happens if only 1 transfusion in group
                    data_sigcolor=data_sigcolor';
                    m = mean(data_sigcolor, 1, 'omitnan')';
                    s = 0;
                else
                    m = mean(data_sigcolor, 1, 'omitnan')';
                    s = std(data_sigcolor, 1, 'omitnan')';
                end
                plot(timeline, m, 'Color', all_colors{icolor}, 'LineWidth', 2, 'DisplayName', displayname);
                %%%%plot(timeline, m, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2,'DisplayName',displayname);
                hold on;
                patch([timeline; timeline(end : -1 : 1)], ...
                    [m + s; m(end : -1 : 1) - s(end : -1 : 1)], all_colors{icolor}, ...
                    'FaceAlpha', alpha, 'EdgeColor', all_colors{icolor},'EdgeAlpha',alpha,'HandleVisibility','off');
                
                plot([0, 0], [min(m - s), max(m + s)], 'k--', 'LineWidth', 2,'HandleVisibility','off');
                
                hold on
                %%% plot([mean(tld.t_stop(idx_include), 'omitnan'), mean(tld.t_stop(idx_include), 'omitnan')], [min(m - s), max(m + s)], 'k--', 'LineWidth', 2);
            
            end
        end
        legend('NumColumns', 3, 'Location', 'northoutside');
    
    %     idx_include = find(sig_quality(isig, :));
    %     m = mean(tld.(var_labels{v})(:, idx_include), 2, 'omitnan');
    %     s = std(tld.(var_labels{v})(:, idx_include), [], 2, 'omitnan');
    %     for ipat = 1:size(signals,3)
    %        thecolor=colors(isig,ipat);
    %        % Signal quality was good and not classified as non-responder
    %        if (tld.sig_quality(isig, ipat) == 1) && all(thecolor ~= [1, 0.41, 0.38])
    %            plot(squeeze(timelines(isig,ipat,:)), ...
    %                 squeeze(signals(isig,ipat,:)), ...
    %                 color=thecolor);
    %            hold on
    %        end
    %     end
    %     % Background color
    %     set(gca,'Color',[1 1 1])
    
        % Looking nice
        set(findall(0, 'type', 'axes'), 'FontName', 'Times', 'Fontsize', 16, 'TickDir', 'out', 'box', 'off', 'linewidth', 2, 'ticklength', [0.01, 0.01])
        orient(fig, 'landscape')
        
        xlabel('Time [hours]')
        ylabel(y_names{isig});
        xlim([timeline(1),timeline(end)])
        set(gcf,'Units','normalized','OuterPosition',[0 0 0.25 0.3]);

        % Save
        exportgraphics(fig, sprintf('%s/%s_vs_time_to_event.jpg', plot_dir, sprintf(var_labels{isig})), 'BackgroundColor', 'None', 'Resolution', 300);
        %print(sprintf('../results/%s_all_patients_%s.jpg',sprintf(var_labels{isig}),dataset_label), '-djpg','-bestfit');
        savefig(fig,sprintf('%s/%s_vs_time_to_event.fig',plot_dir, sprintf(var_labels{isig})) )

    end
    
    %plot piecharts for type of responses
    fig = figure;
    for isig =1:numel(var_labels)
        both=length(find(colors(isig,:)==1));
        increase=length(find(colors(isig,:)==2));
        decrease=length(find(colors(isig,:)==3));
        nochange=length(find(colors(isig,:)==4));
        tot_num=length(find(sig_quality(isig)));
    
        subplot(ceil(numel(var_labels)/3),3,isig)
        
         % Data and labels
         values = [increase, decrease, nochange];
         labels = {'Increase', 'Decrease', 'No Change'};
         percentages = 100 * values / tot_num;
         label_strings = arrayfun(@(n, p) sprintf('%s\n%d (%.1f%%)', labels{n}, values(n), percentages(n)), 1:3, percentages, 'UniformOutput', false);

         subplot(ceil(numel(var_labels)/3), 3, isig)

         % Plot pie
         p = pie(values, label_strings);

         % Set colors
         slice_colors = {
         [0.99, 150/255, 0],   % Increase (orange)
         [0.65, 0.77, 0.90],   % Decrease (blue)
         [1, 0.41, 0.38]       % No change (red)
         };

         for i = 1:2:length(p)             
             set(p(i), 'FaceColor', slice_colors{(i + 1)/2});
         end
        title(y_names{isig});

    end
    set(findall(0, 'type', 'axes'), 'FontName', 'Times', 'Fontsize', 24, 'TickDir', 'out', 'box', 'off', 'linewidth', 2, 'ticklength', [0.01, 0.01])

    set(gcf,'Units','normalized','OuterPosition',[0 0 0.5 0.3]);


    exportgraphics(gcf, sprintf('%s/piechart_all_variables.jpg', plot_dir), 'BackgroundColor', 'none','Resolution',300);

    % save
    print(sprintf('%s/piechart_all_variables.pdf', plot_dir),'-dpdf','-bestfit')
    

%     %plot barcharts for type of responses
fig = figure;

for isig = 1:numel(var_labels)
    % Counts
    both = sum(colors(isig,:) == 1);  % Ignored
    increase = sum(colors(isig,:) == 2);
    decrease = sum(colors(isig,:) == 3);
    nochange = sum(colors(isig,:) == 4);
    tot_num = sum(sig_quality(isig,:));

    % Values and percentages
    values = [increase, nochange, decrease];
    percentages = 100 * values / tot_num;

    % Plot stacked bar
    subplot(ceil(numel(varp_labels)/3), 3, isig)
    b = bar(1, percentages, 'stacked');

    % Set consistent colors
    b(1).FaceColor = [0.99, 150/255, 0];   % Increase (orange)
    b(2).FaceColor = [1, 0.41, 0.38];      % No change (red)
    b(3).FaceColor = [0.65, 0.77, 0.90];   % Decrease (blue)
   

    ylim([0 100])
    ylabel('% of patients')
    title(y_names{isig});
    set(gca, 'XTickLabel', '', 'XTick', []);

    % Add text labels on bars
    y_offset = 0;
    for i = 1:3
        text(1, y_offset + percentages(i)/2, ...
            sprintf('%d\n(%.1f%%)', values(i), percentages(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 12, 'Color', 'k');
        y_offset = y_offset + percentages(i);
    end
end

% Aesthetic settings
set(findall(0, 'type', 'axes'), 'FontName', 'Times', 'Fontsize', 24, ...
    'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'TickLength', [0.01, 0.01]);

set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 0.5 0.5]);

% Export
exportgraphics(gcf, sprintf('%s/barplot_all_patients_%s.jpg', plot_dir, dataset_label), ...
    'BackgroundColor', 'none', 'Resolution', 300);

print(sprintf('%s/barplot_all_variables.pdf', plot_dir), '-dpdf', '-bestfit');

end

end

function idx_full_hb = filter_full_hb(data)
    % Now they always all have the Hb data because we only care about the
    % pre-transfusion
    
    %assert(~sum(isnan(data(1,:)) ~= 0));
    
    idx_full_hb = true(1,size(data,2));% sum(isnan(data),1)==0;
end
% _ EOF____________________________________________________________________
