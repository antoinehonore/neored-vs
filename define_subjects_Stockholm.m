% Defines the Oxford participants from the transfusion project that should 
% be analysed. Note that the base_folder should likely be changed.
%
%
% See also MOVING_AVERAGE, TIME_LOCK_DATA, PLOT_VITAL_SIGNS, CHECK_SIGNALS,
%          FIND_RESPONDERS
% ________________________________________________________________________
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html
%
%                                           (c) Coen Zandvoort, 2024
% ________________________________________________________________________


% get info on infants and recordings
base_folder = '../../matlabdata8/*';
%base_folder = '../data/';
sub_ids = dir(base_folder);


files = cell(1, 1);
file_counter = 1;
for sub = 1 : numel(sub_ids)
    
    ses_ids = dir(fullfile(sub_ids(sub).folder, sub_ids(sub).name, 'RBC'));
    
    for ses = 1 : numel(ses_ids)
        
        file_ids = dir(fullfile(ses_ids(ses).folder, ses_ids(ses).name, '*.mat'));
        
        for f = 1 : numel(file_ids)
            vs = load(fullfile(file_ids(f).folder, file_ids(f).name));
            rates_data = get_rates_data(vs.evt_rates_dose);
            if (vs.birthga/7 <= 32) %&& (vs.evt_start_days < evt_start_days_upper) && (vs.evt_start_days >= evt_start_days_lower)%&& ((vs.evt_sum_dose / vs.weight_kg_pre)>=8) && ((rates_data / vs.weight_kg_pre)>=3)
                files(file_counter, 1) = {fullfile(file_ids(f).folder, file_ids(f).name)};
                %             files(file_counter, 2) = {fullfile(file_ids(f).folder, [file_ids(f).name(1 : strfind(file_ids(f).name, '_vital_signs_data.mat')), 'ibi.mat'])};
                file_counter = file_counter + 1;
            end
        end

    end

end
