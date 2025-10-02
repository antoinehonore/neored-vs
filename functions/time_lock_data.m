function tld = time_lock_data(vs, tld, t_limits, t_baseline, fs)
% - function tld = time_lock_data(vs, tld, t_baseline, fs)
% --time-locks data to transfusion start
%
% Input
% - vs              - Recording-specific vital signs data as struct
% - tld             - Time-locked heart rate, saturation, and respiratory
%                     rate data
% - t_baseline      - Time limits in hours used for baseline correction
%                     (i.e., mean subtraction) defined as double variable.
% - fs              - sampling rate
%
%
% Output
% - tld             - Time-locked heart rate, saturation, and respiratory
%                     rate data as struct.
%
%
% See also PLOT_VITAL_SIGNS, MOVING_AVERAGE, CHECK_SIGNALS, FIND_RESPONDERS
% ________________________________________________________________________
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html
%
%                                           (c) Coen Zandvoort, 2024
% ________________________________________________________________________

% find transfusion start events


idx_events = find(contains(vs.events_vital_signs, 'transfusion') & ...
    contains(vs.events_vital_signs, 'start'));

for e = 1 : numel(idx_events)

    % vital signs data time locking
    [~, idx_data] = min(abs(vs.time - vs.time_events_vital_signs(idx_events(e))));
    time_locked = vs.time - vs.time(idx_data);
    
    %delete nans in time
    index_to_delete=find(isnan(time_locked));
    vs.HR(index_to_delete)=[];
    %vs.fiO2(index_to_delete)=[];
    %vs.Hb(index_to_delete)=[];
    vs.RR(index_to_delete)=[];
    vs.pulse(index_to_delete)=[];
    vs.sats(index_to_delete)=[];
    vs.time(index_to_delete)=[];
    time_locked(index_to_delete)=[];

    %check if any points are sampled twice and remove them
    index_to_delete=find(diff(time_locked)<(fs/3600/2));
    index_to_delete=index_to_delete+1;
    
    vs.HR(index_to_delete)=[];
    vs.RR(index_to_delete)=[];
    vs.pulse(index_to_delete)=[];
    vs.sats(index_to_delete)=[];
    vs.time(index_to_delete)=[];
    time_locked(index_to_delete)=[];
    
    %fs = 0.9766; % Hz (note this is not true but it's a very good approximation). 
    
    % d=diff(vs.time);
    % fs=1/(median(d)*3600);
    tld.time_ref = t_limits(1) : (1 / fs) / 3600  : t_limits(2); % hours
    idx_time = time_locked >= tld.time_ref(1) & time_locked <= tld.time_ref(end);
    idx_first_sample = find(idx_time, 1, 'first');
    [~, idx_start] = min(abs(tld.time_ref - time_locked(idx_first_sample)));
    idx_shift = idx_start : idx_start + sum(idx_time) - 1;
    
    % extract time courses and subtract mean
    idx_baseline = time_locked > t_baseline(1) & time_locked < t_baseline(2);
       
    if strcmp(vs.dataset_label, 'Oxford')
        info = readtable('demographics.xlsx');
        
        for o = 1 : numel(info.RECD_ResearchProject)
            i = contains(vs.path, strcat(info.RECD_ResearchProject(o), info.RECD_ParticipantNumber(o), info.TESTO_TestOccasionNumber(o)), 'IgnoreCase', true);
            if i
                idx = o;
            end
        end
    
        % get demographic data
        tld.BIRTH_Gender{end+1} = info.BIRTH_Gender{idx};
        tld.hb_pre(end+1) = info.StartHB(idx);
        tld.hb_post(end+1) = info.EndHB(idx);
%         tld.fiO2_pre(end+1) = vs.fiO2_pre;
%         tld.fiO2_post(end+1) = vs.fiO2_post;
        tld.PMA(end+1) = info.GA(idx);
        tld.evt_start_pna_days(end+1) = info.PMA(idx) - info.GA(idx);
        tld.evt_start_pma_days(end+1) = info.PMA(idx);
        tld.transfusion_volume(end+1) = info.Volume_kg(idx);
%         tld.transfusion_rate(end+1) = vs.rate;
    
        tld.BIRTH_weight(end+1) = info.BIRTH_BirthWeight(idx);
        tld.TESTOD_CurrentVentilation(end+1) = info.TESTOD_CurrentVentilation(idx);
    
        tld.TESTOD_MostRecentWeight(end+1) = info.TESTOD_MostRecentWeight(idx);
    end

    if strcmp(vs.dataset_label, 'Stockholm')
        tld.BIRTH_Gender{end+1} = vs.sex;
        tld.hb_pre(end+1) = vs.hb_pre;
        tld.hb_post(end+1) = vs.hb_post;
        tld.fiO2_pre(end+1) = vs.fiO2_pre;
        tld.fiO2_post(end+1) = vs.fiO2_post;
        tld.PMA(end+1) = vs.birthga;
        tld.evt_start_pna_days(end+1) = vs.evt_start/24;
        tld.evt_start_pma_days(end+1) = vs.evt_start/24 + vs.birthga;
        tld.transfusion_volume(end+1) = vs.evt_sum_dose;
        tld.transfusion_rate(end+1) = vs.rate;
    
        tld.BIRTH_weight(end+1) = vs.birthweight;
        tld.TESTOD_CurrentVentilation{end+1} = vs.respirator{1,1};
    
        tld.TESTOD_MostRecentWeight(end+1) = vs.weight_kg_pre;
        if vs.weight_kg_pre > 5
            vs.weight_kg_pre
        end
    
    end
    if strcmp(vs.dataset_label, 'Berlin')
        tld.BIRTH_Gender{end+1} = vs.sex;
        tld.hb_pre(end+1) = vs.hb_pre;
        tld.hb_post(end+1) = vs.hb_post;
        tld.fiO2_pre(end+1) = vs.fiO2_pre;
        tld.fiO2_post(end+1) = vs.fiO2_post;
        tld.PMA(end+1) = vs.birthga;
        tld.evt_start_pna_days(end+1) = vs.evt_start/24;
        tld.evt_start_pma_days(end+1) = vs.evt_start/24 + vs.birthga;
        tld.transfusion_volume(end+1) = vs.evt_sum_dose;
        tld.transfusion_rate(end+1) = vs.rate;
    
        tld.BIRTH_weight(end+1) = vs.birthweight;
        tld.TESTOD_CurrentVentilation{end+1} = vs.respirator{1,1};
    
        tld.TESTOD_MostRecentWeight(end+1) = vs.weight_kg_pre;
        if vs.weight_kg_pre > 5
            vs.weight_kg_pre
        end
        
    end
    if strcmp(vs.dataset_label, 'Imperial')
        tld.BIRTH_Gender{end+1} = vs.sex;
        tld.hb_pre(end+1) = vs.hb_pre;
        tld.hb_post(end+1) = vs.hb_post;
         tld.fiO2_pre(end+1) = NaN;
         tld.fiO2_post(end+1) = NaN;
        tld.PMA(end+1) = vs.PMA;
        tld.evt_start_pna_days(end+1) = vs.pna;
        tld.evt_start_pma_days(end+1) = vs.PMA_at_start;
        tld.transfusion_volume(end+1) = vs.evt_sum_dose;
        tld.transfusion_rate(end+1) = vs.rate;
        % 
        tld.BIRTH_weight(end+1) = vs.birthweight;
         tld.TESTOD_CurrentVentilation{end+1} = NaN; %vs.respirator{1,1};
        % 
         tld.TESTOD_MostRecentWeight(end+1) = NaN; %vs.weight_kg_pre;
        % if vs.weight_kg_pre > 5
        %     vs.weight_kg_pre
        % end
        % 
    end
    tmp = vs.HR(idx_time);
    if ~strcmp(vs.dataset_label,'Stockholm')
        vs.HR(vs.HR == 0) = NaN;
        tld.hr(:, tld.counter) = NaN(size(tld.time_ref, 2), 1);
        tmp = interp1(time_locked(idx_time),vs.HR(idx_time),tld.time_ref);
    end
    
    if strcmp(vs.dataset_label, 'Stockholm')
        tld.hr(idx_shift, tld.counter) = tmp - mean(vs.HR(idx_baseline), 'omitnan'); %interpolate to account for sometimes sporadic sampling rate (particularly in Berlin data)
    else
        tld.hr(:, tld.counter) = tmp - mean(vs.HR(idx_baseline), 'omitnan'); %interpolate to account for sometimes sporadic sampling rate (particularly in Berlin data)
    end
    
    tld.hr(tld.hr(:, tld.counter) == 0, tld.counter) = NaN;
    
    tmp = vs.sats(idx_time);
    if ~strcmp(vs.dataset_label,'Stockholm')
        vs.sats(vs.sats == 0) = NaN;
        tld.sats(:, tld.counter) = NaN(size(tld.time_ref, 2), 1);
        tmp = interp1(time_locked(idx_time),vs.sats(idx_time),tld.time_ref);
    end
    if strcmp(vs.dataset_label, 'Stockholm')
        tld.sats(idx_shift, tld.counter) = tmp - mean(vs.sats(idx_baseline), 'omitnan');
    else
        tld.sats(:, tld.counter) = tmp - mean(vs.sats(idx_baseline), 'omitnan');
    end

    tld.sats(tld.sats(:, tld.counter) == 0, tld.counter) = NaN;

    tmp = vs.RR(idx_time);
    if ~strcmp(vs.dataset_label,'Stockholm')
        vs.RR(vs.RR == 0) = NaN;
        tld.rr(:, tld.counter) = NaN(size(tld.time_ref, 2), 1);
        tmp = interp1(time_locked(idx_time),vs.RR(idx_time),tld.time_ref) ;
    end
    if strcmp(vs.dataset_label, 'Stockholm')
        tld.rr(idx_shift, tld.counter) = tmp - mean(vs.RR(idx_baseline), 'omitnan');
    else
        tld.rr(:, tld.counter) = tmp - mean(vs.RR(idx_baseline), 'omitnan');
    end
    
    tld.rr(tld.rr(:, tld.counter) == 0, tld.counter) = NaN;
    
    %tmp = vs.fiO2(idx_time);
    %if ~strcmp(vs.dataset_label,'Stockholm')
    %    vs.fiO2(vs.fiO2 == 0) = NaN;
    %    tld.fiO2(:, tld.counter) = NaN(size(tld.time_ref, 2), 1);
    %    tmp = interp1(time_locked(idx_time),vs.fiO2(idx_time),tld.time_ref) ;
    %end
    %tld.fiO2(idx_shift, tld.counter) = tmp - mean(vs.fiO2(idx_baseline), 'omitnan');
    %tld.fiO2(tld.fiO2(:, tld.counter) == 0, tld.counter) = NaN;
    
    %tmp = vs.Hb(idx_time);
    %if ~strcmp(vs.dataset_label,'Stockholm')
    %    vs.Hb(vs.Hb == 0) = NaN;
    %    tld.Hb(:, tld.counter) = NaN(size(tld.time_ref, 2), 1);
    %    tmp = interp1(time_locked(idx_time),vs.Hb(idx_time),tld.time_ref) ;
    %end
    %tld.Hb(idx_shift, tld.counter) = tmp; %%%- mean(vs.Hb(idx_baseline), 'omitnan');
    %tld.Hb(tld.Hb(:, tld.counter) == 0, tld.counter) = NaN;
    

    % time-lock the inter-breath intervals
    if isfield(vs, 'ibi')
        vs.ibi.ibi_time_locked = vs.ibi.ibi_time ./ 3600;
        [~, idx_data] = min(abs(vs.ibi.ibi_time_locked - vs.time_events_vital_signs(idx_events(e))));
        vs.ibi.ibi_time_locked = vs.ibi.ibi_time_locked - vs.ibi.ibi_time_locked(idx_data);
        idx_time = vs.ibi.ibi_time_locked >= tld.time_ref(1) & vs.ibi.ibi_time_locked <= tld.time_ref(end);
        
        tld.ibi_locked{tld.counter} = vs.ibi.ibi(idx_time);
        tld.ibi_time_locked{tld.counter} = vs.ibi.ibi_time_locked(idx_time);
    end


    % find if transfusion stop is available
    idx_event_stop = find(contains(vs.events_vital_signs(idx_events(e) : end), 'transfusion') & contains(vs.events_vital_signs(idx_events(e) : end), 'stop')) + idx_events(e) - 1;
    if ~isempty(idx_event_stop)
        if length(idx_event_stop)==length(idx_events)
            tld.t_stop(tld.counter, 1) = vs.time_events_vital_signs(idx_event_stop(e)) - vs.time_events_vital_signs(idx_events(e));
        else
            tld.t_stop(tld.counter, 1) = NaN; %for now set to NaN if not enough end markers
        end
    else
        tld.t_stop(tld.counter, 1) = NaN;
    end
   
   
    tld.counter = tld.counter + 1;

end

end

% _ EOF____________________________________________________________________