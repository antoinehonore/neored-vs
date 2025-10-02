function [T, clusters_start_original, clusters_end_original, pvals_mass, clusters_size_original, max_clusters_size] = montecarlo_fornewdata_differentsizes(dataPre, dataPost, fs, presec, postsec)

%input - dataPre - baseline data
%       dataPost - stimulus-evoked data (same order as dataPre), for this
%       function the length of the pre and post data can be different
%       lengths
%       fs - sampling rate
%       presec - prestimulus epoch length in seconds
%       postsec - poststimulus epoch length in seconds

% Caroline Hartley


% For reproducibility
rng(1)


% Instead of taking differences use the whole pre data as one distribution,
% compare to each point in the post
% take the median t statistic

n1=size(dataPost,1); %trials
n2=size(dataPost,2); %time points
tstore=zeros(n1,n2);
for timepoint=1:n2
for trial=1:n1
    if ~any(~isnan(dataPre(trial,:))) && isnan(dataPost(trial,timepoint)); continue; end
    [~,~,~,stats] = ttest(dataPre(trial,:),dataPost(trial,timepoint));
    tstore(trial,timepoint) = stats.tstat;
end
end

allt=median(tstore,"omitnan");

% threshold: degrees of freedom = num subjects - 1
thres = tinv(0.975,(n1-1));
[clusters_no,clusters_start_original,clusters_end_original,clusters_length_original] = seq_cluster(abs(allt),thres,1,presec,postsec,fs,0);

clusters_size_original = zeros(clusters_no,1);
for i = 1:clusters_no
    clusters_size_original(i) = sum(abs(allt(clusters_start_original(i):clusters_end_original(i))));
end


% only if a cluster is found, compare with random sample
if ~isempty(clusters_start_original)
    
    %join pre and post data so can permute all of it together
    dataall=[dataPre,dataPost];
    s1=size(dataPre,2);

    perm_no = 1000;
    max_clusters=zeros(perm_no,1);
    max_clusters_size=zeros(perm_no,1);

    % generate permutations
    Pset = palm_quickperms(ones(size(dataall,1),1), [], perm_no, 0, 1); % perms
    Pset = sign(Pset); % signs
    
    for j = 1:perm_no % first permutation is the true data
     %   display([num2str(j), ' out of ', num2str(perm_no)])
        
        %  A) Grab flips
        signs = Pset(:,j);
        %  B) Multiply differences by this
        permData = dataall .* signs; 
        %  C) Split data and Calculate paired t-statistic
        dataPre_shuff=permData(:,1:s1);
        dataPost_shuff=permData(:,s1+1:end);

        tstore=zeros(n1,n2);
        for timepoint=1:n2
        for trial=1:n1
            if ~any(~isnan(dataPre_shuff(trial,:))) && isnan(dataPost_shuff(trial,timepoint)); continue; end
            [~,~,~,stats] = ttest(dataPre_shuff(trial,:),dataPost_shuff(trial,timepoint));
            tstore(trial,timepoint) = stats.tstat;
        end
        end

        allt_perm=median(tstore,"omitnan");

        %  D) Find clusters
        [clusters_no,clusters_start,clusters_end,clusters_length] = seq_cluster(abs(allt_perm),thres,1,presec,postsec,fs,0);
        
        if ~isempty(clusters_length)
            max_clusters(j) = max(clusters_length); % max cluster length for this permutation
            clusters_size = zeros(clusters_no,1);
            for c = 1:clusters_no
                clusters_size(c) = sum(abs(allt_perm(clusters_start(c):clusters_end(c)))); % cluster mass (sum of t-stats)
            end
            max_clusters_size(j) = max(clusters_size); % max cluster mass for this permutation
        end
    end

    % calculate p-values: for each cluster in the true data, count how many
    % clusters in the permutations were of the same or larger mass
    pvals_mass = NaN(length(clusters_start_original),1);
    for i = 1:length(clusters_length_original)
        pvals_mass(i) = sum(max_clusters_size >= clusters_size_original(i)) /perm_no;
    end

    % Let's add a funky histogram to view this
%     numbins = floor(perm_no / 2); 
%     figure; 
%     for i = 1:length(clusters_size_original)
%         subplot(2,1,1); hold on; 
%         histogram(max_clusters_size,numbins); 
%         plot([clusters_size_original(i),clusters_size_original(i)],[0 0.1*perm_no])
%         title('Cluster size')
%     end

else
    pvals_mass = NaN;
    max_clusters_size = NaN;
    clusters_start_original = NaN; 
    clusters_end_original = NaN;
    clusters_size_original = NaN;
end
T = table((clusters_start_original/fs), (clusters_end_original/fs), pvals_mass, 'VariableNames', {'Cluster_start','Cluster_end','p_value'});
%T = table((clusters_start_original/fs),(clusters_end_original/fs),pvals_mass);
%T.Properties.VariableNames={'Cluster_start','Cluster_end','p_value'};
%T
disp(T)