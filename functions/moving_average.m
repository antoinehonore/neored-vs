function tld = moving_average(tld, t_window, t_overlap)
% - function tld = moving_average(tld, t_window, t_overlap) -- computes a
% moving average
%
%
% Input
% - tld             - Time-locked heart rate, saturation, and respiratory
%                     rate data
% - t_window        - Window size in hours defined as double variable.
% - t_overlap       - Overlap between consecutive windows in hours defined
%                     as double.
%
%
% Output
% - tld             - Time-locked heart rate, saturation, and respiratory
%                     rate data 
%
%
% See also PLOT_VITAL_SIGNS, TIME_LOCK_DATA, CHECK_SIGNALS, FIND_RESPONDERS
% ________________________________________________________________________
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html
%
%                                           (c) Coen Zandvoort, Caroline Hartley 2024
% ________________________________________________________________________


% compute moving averages
n_windows = round((diff(tld.time_ref([1, end])) - t_window + t_overlap) / t_overlap);
var_names = {'hr', 'sats', 'rr'};

for f=1:tld.counter-1
for v = 1 : numel(var_names)

    for w = 1 : n_windows

        % get indexes
        start_stop = [0, t_window] + tld.time_ref(1) + t_overlap * (w - 1);
        idx_time = tld.time_ref > start_stop(1) & tld.time_ref <= start_stop(2);

        tld.(strcat(var_names{v}, '_mean'))(w, f) = mean(tld.(var_names{v})(idx_time, f), 'omitnan');
        tld.(strcat(var_names{v}, '_std'))(w, f) = std(tld.(var_names{v})(idx_time, f), 'omitnan');
        tld.t_average(w, 1) = start_stop(end);

        if isfield(tld, 'ibi_locked')
            
            idx_time = tld.ibi_time_locked > start_stop(1) & tld.ibi_time_locked <= start_stop(2);
            if sum(idx_time) ~= 0
                tld.resp_rate(w, f) = 60 ./ mean(tld.ibi_locked(idx_time), 'omitnan'); % breaths/min
                tld.ibi_std(w, f) = std(tld.ibi_locked(idx_time), 'omitnan');
                tld.apnoea_rate_5_sec(w, f) = sum(tld.ibi_locked(idx_time) > 5); % apnoeas/hour
                tld.apnoea_rate_10_sec(w, f) = sum(tld.ibi_locked(idx_time) > 10); % apnoeas/hour
            else
                tld.resp_rate(w, f) = NaN;
                tld.ibi_std(w, f) = NaN;
                tld.apnoea_rate_5_sec(w, f) = NaN;
                tld.apnoea_rate_10_sec(w, f) = NaN;
            end

        end

    end

end


end

end

% _ EOF____________________________________________________________________