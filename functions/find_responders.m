function [signals, colors,all_colors,signal_increase,signal_decrease,all_colors_names] = find_responders(tld, std_threshold, sample_threshold,plot_subject_specific)
% - function find_responders(tld, std_threshold, sample_threshold) -- 
% finds responders using pre-defined threshold
%
%
% Input
% - tld                     - Time-locked heart rate, saturation, and 
%                             respiratory rate data
% - std_threshold           - Threshold used to define (non-)responders;
%                             it should be exceeded post-transfusion start 
%                             to be classified as responder
% - sample_threshold        - Threshold used to define (non-)responders;
%                             post-transfusion response should exceed this
%                             number of consecutive samples to be 
%                             classified as responder.
%
%
% See also MOVING_AVERAGE, TIME_LOCK_DATA, PLOT_VITAL_SIGNS, CHECK_SIGNALS
% ________________________________________________________________________
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html
%
%                                           (c) Coen Zandvoort, Antoine Honore, Caroline Hartley 2024
% ________________________________________________________________________
if nargin < 4; sample_threshold = 6; plot_subject_specific=1; end

% find responders
var_labels = {'hr_mean', 'sats_mean', 'rr_mean', 'hr_std', 'sats_std', 'rr_std'};
if isfield(tld, 'ibi_locked'); var_labels = [var_labels, {'resp_rate', 'ibi_std', 'apnoea_rate_5_sec', 'apnoea_rate_10_sec'}]; end
y_names = {'Mean HR [bpm]', 'Mean Sats [%]', 'Mean RR [breaths/min]', 'Std HR [bpm]', 'Std Sats [%]', 'Std RR [breaths/min]'};
if isfield(tld, 'ibi_locked'); y_names = [y_names, {'Resp rate [breaths/min]', 'Std IBI [sec]', 'Apnoea rate 5 sec [times/hour]', 'Apnoea rate 10 sec [times/hour]'}]; end

close all
enable_plot = false;
T = size(tld.t_average,1);
n_evts = tld.counter-1;

colors = zeros(numel(var_labels), n_evts);
signals = zeros(numel(var_labels), n_evts, T);
signals(:)=NaN;
signal_increase = zeros(numel(var_labels), n_evts);
signal_decrease = zeros(numel(var_labels), n_evts);


% all_colors = {'green','red','blue','black'};
%all_colors = {'green',[0.99, 0.99, 0.59],[0.65, 0.77, 0.90],[1, 0.41, 0.38]};
all_colors = {'green',[0.99, 150/255, 0],[0.65, 0.77, 0.90],[1, 0.41, 0.38]};
all_colors_names = {'Both', 'Increase', 'Decrease', 'Stable'};


for v = 1 : numel(var_labels)
    % find indexes pre- and post-transfusion start
    idx_pre = tld.t_average <= -0.5;
    std_pre = std(tld.(var_labels{v})(idx_pre, :), 'omitnan');
    idx_post = tld.t_average >= 0;

    for r = 1 : numel(std_pre)

        % check signal quality
        if tld.sig_quality(v, r) == 0; colors(v,r)=NaN; signal_increase(v,r)=NaN; signal_decrease(v,r)=NaN; continue; end

        % find if the post-transfusion signal exceeded the std threshold
        idx_exceed = (find(abs(tld.(var_labels{v})(idx_post, r)) > std_pre(r) .* std_threshold))';
        idx_consecutive = diff(idx_exceed) == 1;
        %check that no consecutive points are a switch
        p=find(idx_post);
        s=sign(tld.(var_labels{v})(p(idx_exceed), r));
        switch_sign=find(diff(s));
        idx_consecutive(switch_sign)=0;
        idx_start_stop = find([false, idx_consecutive] ~= [idx_consecutive, false]);
        gaps = find(idx_start_stop(2 : 2 : end) - idx_start_stop(1 : 2 : end - 1) >= sample_threshold);
        
        %calculate amount of signal above and below the threshold
        cutsignalabove = tld.(var_labels{v})(idx_post, r) - (std_pre(r) .* std_threshold);  %cut signal so get sum above threshold (i.e. subtracting variability in baseline)
        signal_increase(v,r) = sum(cutsignalabove(find(cutsignalabove>0)))/length(cutsignalabove);
        
        cutsignalbelow = tld.(var_labels{v})(idx_post, r) - (-std_pre(r) .* std_threshold); 
        signal_decrease(v,r) = sum(abs(cutsignalbelow(find(cutsignalbelow<0))))/length(cutsignalbelow);

        
        % change color 0 - which identifies if a baby is a responder or not
        idx_start = find(idx_post, 1, 'first') + idx_exceed(idx_start_stop(2 * gaps - 1));

        if any(tld.(var_labels{v})(idx_start, r) < 0) && any(tld.(var_labels{v})(idx_start, r) > 0)
            color = 1;%'green';%[0.47, 0.87, 0.47];
        elseif any(tld.(var_labels{v})(idx_start, r) > 0)
            %color = [0.99, 0.99, 0.59];%Yellow
            color = 2;%'black';%[0,0,0]; %Black
        elseif any(tld.(var_labels{v})(idx_start, r) < 0)
            color = 3;%'blue';%[0.65, 0.77, 0.90];%blue
        else
            color = 4;%'red';%[1, 0.41, 0.38];%red
        end
        
        % Save data in tensors that are returned
        colors(v,r) = color;
        signals(v,r,:) = tld.(var_labels{v})(:, r);
        
        % % plot subject-specific responses
        if plot_subject_specific
            figure(v * 10); %('Visible','off')
            po = get(gcf, 'position');
            set(gcf, 'position', [po(1:2), 1200, 900], 'name', sprintf('recording_responses_%s', var_labels{v}));
            
            % plot mean and standard deviations
            s1 = numel(std_pre);
            s2 = ceil(sqrt(s1));
            subplot(s2, ceil(s1 / s2), r);
            
            y_lim = max(abs(min(tld.(var_labels{v})(:))), max(tld.(var_labels{v})(:)));
            
            plot(tld.t_average, tld.(var_labels{v})(:, r), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
            hold on;
            plot([0, 0], [-y_lim, y_lim], 'k--', 'LineWidth', 2);
            plot([tld.t_stop(r), tld.t_stop(r)], [-y_lim, y_lim], 'k--', 'LineWidth', 1.5);
            plot([tld.t_average(1), tld.t_average(end)], [0, 0], 'k', 'LineWidth', 1.5)
            plot([tld.t_average(1), tld.t_average(end)], [std_pre(r) .* std_threshold, std_pre(r) .* std_threshold], 'k--', 'LineWidth', 1.5)
            plot([tld.t_average(1), tld.t_average(end)], -[std_pre(r) .* std_threshold, std_pre(r) .* std_threshold], 'k--', 'LineWidth', 1.5)
            
            xlabel('Time [hours]')
            ylabel(y_names{v});
            xlim([tld.t_average([1, end])])
            ylim([-y_lim, y_lim])
    
            % change background color
            idx_start = find(idx_post, 1, 'first') + idx_exceed(idx_start_stop(2 * gaps - 1));
            
            if any(tld.(var_labels{v})(idx_start, r) < 0) && any(tld.(var_labels{v})(idx_start, r) > 0)
                set(gca, 'Color', [0.47, 0.87, 0.47], 'GridAlpha', 0.3); 
            elseif any(tld.(var_labels{v})(idx_start, r) > 0)
                set(gca, 'Color', [0.99, 0.99, 0.59], 'GridAlpha', 0.3); 
            elseif any(tld.(var_labels{v})(idx_start, r) < 0)
                set(gca, 'Color', [0.65, 0.77, 0.90], 'GridAlpha', 0.3); 
            else
                set(gca, 'Color', [1, 0.41, 0.38], 'GridAlpha', 0.3); 
            end
            % 
            % print(sprintf('../matlabplot/responders/patient_%d_%s',r,var_labels{v}),'-dpdf','-bestfit')
            % close all
        end
    end
end

end

% _ EOF____________________________________________________________________