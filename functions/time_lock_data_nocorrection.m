function tld = time_lock_data_nocorrection(vs, tld, t_limits, fs)
% --time-locks data to transfusion start
% NB use time_lock_data fundction if you want to baseline correct
%
% Input
% - vs              - Recording-specific vital signs data as struct
% - tld             - Time-locked heart rate, saturation, and respiratory
%                     rate data
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
% Edited by Caroline Hartley
% ________________________________________________________________________

% find transfusion start events
idx_events = find(contains(vs.events_vital_signs, 'transfusion') & ...
    contains(vs.events_vital_signs, 'start'));
max_gap_length = 15;

for e = 1 : numel(idx_events)

    % vital signs data time locking
    [~, idx_data] = min(abs(vs.time - vs.time_events_vital_signs(idx_events(e))));
    time_locked = vs.time - vs.time(idx_data);
    
    
    %fs = 0.9766; % Hz (note this is not true but it's a very good approximation). 
    tld.time_ref = t_limits(1) : (1 / fs) / 3600  : t_limits(2); % hours
    idx_time = time_locked >= tld.time_ref(1) & time_locked <= tld.time_ref(end);
    idx_first_sample = find(idx_time, 1, 'first');
    [~, idx_start] = min(abs(tld.time_ref - time_locked(idx_first_sample)));
    idx_shift = idx_start : idx_start + sum(idx_time) - 1;


    % extract time courses

    vs.HR(vs.HR == 0) = NaN;
    tld.hr(:, tld.counter) = NaN(size(tld.time_ref, 2), 1);
    tld.hr(idx_shift, tld.counter) = vs.HR(idx_time);
    tld.hr(tld.hr(:, tld.counter) == 0, tld.counter) = NaN;
    tld.hr(:,tld.counter) = interp_max_interval(tld.hr(:,tld.counter),max_gap_length);

    vs.sats(vs.sats == 0) = NaN;
    tld.sats(:, tld.counter) = NaN(size(tld.time_ref, 2), 1);
    tld.sats(idx_shift, tld.counter) = vs.sats(idx_time);
    tld.sats(tld.sats(:, tld.counter) == 0, tld.counter) = NaN;
    tld.sats(:,tld.counter) = interp_max_interval(tld.sats(:,tld.counter),max_gap_length);

    vs.RR(vs.RR == 0) = NaN;
    tld.rr(:, tld.counter) = NaN(size(tld.time_ref, 2), 1);
    tld.rr(idx_shift, tld.counter) = vs.RR(idx_time);
    tld.rr(tld.rr(:, tld.counter) == 0, tld.counter) = NaN;
    tld.rr(:,tld.counter) = interp_max_interval(tld.rr(:,tld.counter),max_gap_length);

    if isfield(vs, 'ibi')
        vs.ibi.ibi_time_locked = vs.ibi.ibi_time ./ 3600;
        [~, idx_data] = min(abs(vs.ibi.ibi_time_locked - vs.time_events_vital_signs(idx_events(e))));
        vs.ibi.ibi_time_locked = vs.ibi.ibi_time_locked - vs.ibi.ibi_time_locked(idx_data);
        idx_time = vs.ibi.ibi_time_locked >= tld.time_ref(1) & vs.ibi.ibi_time_locked <= tld.time_ref(end);
        
        tld.ibi_locked = vs.ibi.ibi(idx_time);
        tld.ibi_time_locked = vs.ibi.ibi_time_locked(idx_time);
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