% Defines the Imperial participants from the transfusion project that should 
% be analysed. 
% Note that the base_folder should likely be changed.
%
%
% See also MOVING_AVERAGE, TIME_LOCK_DATA, PLOT_VITAL_SIGNS, CHECK_SIGNALS,
%          FIND_RESPONDERS
% ________________________________________________________________________
%
% This file is largely based on the file released under the terms of the 
% GNU General Public License, version 3.
% See http://www.gnu.org/licenses/gpl.html
%
%                                           (c) Coen Zandvoort, 2024
% ________________________________________________________________________
%
% get info on infants and recordings
base_folder = [strrep(pwd, '\','/') '/MATLAB_data_Imperial/'];
if startsWith(base_folder, 'C:')
    base_folder=base_folder(3:end)
end

files = {}; % Initialize as an empty cell array
file_counter = 0; % Start the counter at 0
excluded_file_counter = 0; % Start the counter at 0
skipped_file_counter = 0; % Start the counter at 0

file_ids = dir(fullfile(base_folder, 'Imperial_*'));

excluded_subjects = {};
included_subjects = {};

for f = 1 : numel(file_ids)
    file_path = fullfile(file_ids(f).folder, file_ids(f).name);
    
    % Load .mat file
    vs = load(file_path);
    
    % Basic check for required fields
    if isfield(vs, 'PMA') && isfield(vs, 'pna') && isfield(vs, 'subject')
        % Apply inclusion criteria
        if (vs.PMA / 7 <= 32) && (vs.pna <= evt_start_days_upper) && (vs.pna > evt_start_days_lower)

            file_counter = file_counter + 1;
            files{file_counter, 1} = file_path;
            included_subjects{end + 1} = vs.subject;
        else
            excluded_subjects{end + 1} = vs.subject;
            excluded_file_counter = excluded_file_counter + 1;
        end
    else
        warning('Skipping file (missing PMA, pna, or subject): %s', file_path);
        skipped_file_counter = skipped_file_counter + 1;
    end
end

% Determine unique subject counts
included_unique = numel(unique(included_subjects));
excluded_unique = numel(setdiff(unique(excluded_subjects), unique(included_subjects)));

% Display summary
fprintf('\n--- Inclusion Summary ---\n');
fprintf('Total files processed:      %d\n', numel(file_ids));
fprintf('Files skipped:              %d\n', skipped_file_counter);
fprintf('Files included:             %d\n', file_counter);
fprintf('Files excluded:             %d\n', excluded_file_counter);
fprintf('Unique subjects included:   %d\n', included_unique);
fprintf('Unique subjects excluded:   %d\n\n', excluded_unique);
