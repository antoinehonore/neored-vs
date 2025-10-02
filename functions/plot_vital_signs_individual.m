function plot_vital_signs_individual(files, t_limits, t_baseline,baseline_correct)
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
% - t_baseline      - Time limits in hours used for baseline correction
%                     (i.e., mean subtraction) defined as double variable.
% - baseline_correct - 0 if don't want to baseline correct, 1 if you do
%
%
%
% Example usage
% - define_subjects_Oxford (example script for data structure).
% - t_limits = [-6, 12]; % hours
% - t_baseline = [-6, 0]; % hours
% - plot_vital_signs_individual(files, t_limits,t_baseline)
%
% ________________________________________________________________________
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html
%
%                                           (c) Caroline Hartley, 2024
% ________________________________________________________________________


% check the availability of the functions folder
if exist(fullfile(pwd, 'functions'), 'dir')
    addpath(fullfile(pwd, 'functions'))
else
    error('check if "functions"-folder is part of the current folder')
end


% pre-allocate some struct and fields
tld = struct;
tld.counter = 1;
toremove=[];

var_labels = {'hr', 'sats', 'rr'};
for f = 1 : size(files, 1)

    % load vital signs and inter-breath interval data
    fprintf('reading:\t"%s"\n', files{f, 1})
     newstr=cell2mat(extractBefore(files(f),'_vital_signs_data'));
    toID=newstr(end-8:end);

    vs = load(files{f, 1});
    if size(files, 2) > 1
        if ~isempty(files{f, 2})
            % load inter-breath intervals
            vs.ibi = load(files{f, 2});
        end
    end


    % find transfusion starts - will skip data set if it can't find
    % transfusion start
    if ~isfield(vs, 'events_vital_signs'); disp('no data'); toremove=[toremove;f]; continue; end
    if ~any(find(contains(vs.events_vital_signs, 'transfusion') & ...
        contains(vs.events_vital_signs, 'start'))); disp('no transfusion marker'); toremove=[toremove;f]; continue; end


    % time-lock data
    if baseline_correct
        tld = time_lock_data(vs, tld, t_limits, t_baseline);
    else
        tld = time_lock_data_nocorrection_withtitle(toID, vs, tld, t_limits);
    end
    

end

files(toremove)=[];

% check vital signs signals visually (note that if sig_check is set to
% false all data will be included). 
if isfield(tld, 'ibi_locked'); var_labels = [var_labels, {'resp_rate', 'ibi_std', 'apnoea_rate_5_sec', 'apnoea_rate_10_sec'}]; end


for v = 1 : numel(var_labels)

    x = tld.(sprintf('%s', var_labels{v}));
    time = tld.time_ref;
    events = tld.t_stop;
    figure; plotcount=1;

    for k=1:size(x,2)
        if plotcount==10
            figure;
            plotcount=1;
        end

        subplot(3,3,plotcount)
        plot(time,x(:,k))
        xlim(t_limits)
        hold on
        plot([0,0],ylim,'k--')
        plot([events(k),events(k)],ylim,'k--')
        
        title(tld.toID(k,:))
        ylabel(var_labels{v})
        plotcount=plotcount+1;

    end

end





end

% _ EOF____________________________________________________________________