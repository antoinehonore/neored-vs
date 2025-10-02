% d='/Users/carolinehartley/Desktop/transfusion/';
% files=dir(strcat(d,"*_numerics.csv"));
% 
% for i=1:length(files)
%     files(i).name
%     load_monitor_data(files(i).name, d, d)
% end

%% find data - must run this section first
% define_subjects_Oxford

%% plot raw traces around transfusion
% t_limits = [-6, 12]; % hours
% t_baseline = [-6, 0]; % hours
% plot_vital_signs_individual(files, t_limits, t_baseline,0) %NB change last number to 0 if don't want to baseline correct or 1 if you do
clear all
close all

% Default
dataset_label = 'Oxford';
load_other_centers = true;
cluster_analysis = true;

[ret, name] = system('hostname');
name = strip(name);
if strcmp(name, 'cmm0958')
    addpath('lib/palm/palm-alpha119/')
    dataset_label = 'Stockholm';
end

if contains(name, 'Mac')
    addpath('lib/palm/palm-alpha119/')
    dataset_label = 'Imperial';
    load_other_centers = false;
end

%% plot average and find responders and non-responders

if exist(fullfile(pwd, 'functions'), 'dir')
    addpath(fullfile(pwd, 'functions'))
else
    error('check if "functions"-folder is part of the current folder')
end


t_window = 1; % hours
t_overlap = 0.5; % hours
t_limits = [-12, 12]; % hours
t_baseline = [-6, 0]; % hours
check_sig = true;
std_threshold = 1.0;


if contains(name, 'Mac')
    preprocessed_data_folder = 'cached';
    addpath(fullfile(pwd, 'cached'))
else
    preprocessed_data_folder = '../preprocessed';
end

if ~exist(preprocessed_data_folder, 'dir')
    error('check if "preprocessed"-folder is part of the current folder')
end

sig_quality_fname = sprintf('%s/sig_quality_%s.mat', preprocessed_data_folder, dataset_label);
transfusion_rawdata_fname = strcat(preprocessed_data_folder,'/rawtransfusiondata.mat');

plot_dir_stem = 'results_Stockholm';
if load_other_centers
    plot_dir_stem=strcat(plot_dir_stem,'_Oxford_Imperial');
end

%%  Run analysis with several cutoffs on the age at event start.
for ianalysis = 3:3
    % Remove cache 
    delete(sig_quality_fname)
    delete(transfusion_rawdata_fname)
    
    % Choice of cutoff & output directory
    if ianalysis==1 % all events
        evt_start_days_upper = 1000000000;
        evt_start_days_lower = 0;
        plot_dir = strcat('../',plot_dir_stem);%results_Stockholm_Oxford_Imperial';
        
    elseif ianalysis==2 % events above 2w only
        evt_start_days_upper = 1000000000;
        evt_start_days_lower = 2*7;
        plot_dir = strcat('../',plot_dir_stem,'_above_2w');
    
    elseif ianalysis==3 % events below 2w only
        evt_start_days_upper = 2*7;
        evt_start_days_lower = 0;
        plot_dir = strcat('../',plot_dir_stem,'_below_2w');
    end
    
    mkdir(plot_dir)

    %% Find the relevant data files
    if strcmp(dataset_label, 'Stockholm')
        % Uses the variables evt_start_days_upper and evt_start_days_lower
        define_subjects_Stockholm
    end
    
    if strcmp(dataset_label, 'Oxford')
        % Uses the variables evt_start_days_upper and evt_start_days_lower
        define_subjects_Oxford
    end
    if strcmp(dataset_label, 'Berlin')
        % Uses the variables evt_start_days_upper and evt_start_days_lower
        define_subjects_Berlin
    end
    if strcmp(dataset_label, 'Imperial')
        % Uses the variables evt_start_days_upper and evt_start_days_lower
        define_subjects_Imperial
    end
    
    [signals, responders,signal_increase,signal_decrease,var_labels,tld] = plot_vital_signs(files, ...
        t_limits, t_window, t_overlap, t_baseline, check_sig, dataset_label, std_threshold, plot_dir, ...
        sig_quality_fname, transfusion_rawdata_fname,evt_start_days_lower,evt_start_days_upper,load_other_centers,cluster_analysis);
    
    %% Create table of results and save as excel files
    increase_responder=zeros(size(responders));
    decrease_responder=zeros(size(responders));
    increase_responder(find(responders==2))=1;
    increase_responder(find(responders==1))=1;
    decrease_responder(find(responders==3))=1;
    decrease_responder(find(responders==1))=1;
    
    % Create one table for each variable of interest
    for v=1:length(var_labels)
        T=table(tld.subjectid, tld.studyid, tld.pre_post_hb(1,:)',tld.pre_post_hb(2,:)',tld.TESTOD_CurrentVentilation',tld.TESTOD_MostRecentWeight', tld.PMA',tld.evt_start_pna_days', tld.BIRTH_Gender',increase_responder(v,:)',decrease_responder(v,:)',signal_increase(v,:)',signal_decrease(v,:)', ...
            'VariableNames',{'Subject ID';'Study ID';'StartHB';'EndHB';'TESTOD_CurrentVentilation'; 'TESTOD_MostRecentWeight';'PMA';'PNA';'BIRTH_Gender';'Increase Response';'Decrease Response';'Average Signal Increase';'Average Signal Decrease'});
        excel_fname = strcat(plot_dir,'/tables/', var_labels(v), sprintf('_%s.xlsx',dataset_label));
        mkdir(strcat(plot_dir,'/tables/'))
        writetable(T, cell2mat(excel_fname));
    end
    
    close all
    
    plot_binary_responder
    close all

    %investigate_correlations
    %close all
end



%% Berlin

% define_subjects_Berlin;
% t_window = 1; % hours
% t_overlap = 0.5; % hours
% t_limits = [-6, 12]; % hours
% t_baseline = [-6, 0]; % hours
% check_sig = true;
% dataset_label = 'Berlin';
% std_threshold = 3;
% ibi_boolean = false;
% plot_vital_signs(files, t_limits, t_window, t_overlap, t_baseline, check_sig, dataset_label, std_threshold)