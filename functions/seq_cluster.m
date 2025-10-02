% for a sequence find clusters of values less than a threshold

function [clusters_no,clusters_start,clusters_end,clusters_length]=seq_cluster(series,thres,above_thres,presec,postsec,fs,plot_fig)

% Caroline Hartley May 2016

% output:   clusters_no - number of clusters
%           clusters_start - start value of clusters
%           clusters_end - end time point of clusters
%           clusters_length - number of points in each cluster
% input:    series - time series
%           thres - threshold
%           above_thres=1 or 0 - 1 indicates find points above the
%                   threshold
%           plot_fig = 1 or 0


if above_thres
    below_thres=find(series>thres);
else
   below_thres=find(series<thres);
end

if ~isempty(below_thres)
    clusters_start=below_thres(1);
    clusters_end=[];

    for i=2:length(below_thres)

        if floor(below_thres(i)-below_thres(i-1))>1
            clusters_start=[clusters_start;below_thres(i)];
            clusters_end=[clusters_end;below_thres(i-1)];
        end
    end

    clusters_end=[clusters_end; below_thres(end)];

    clusters_length=clusters_end-clusters_start+1;
    clusters_no=length(clusters_start);
else
    clusters_start=[];
    clusters_end=[];
    clusters_no=[];
    clusters_length=[];
end

if plot_fig
    time = presec:1/fs:postsec;
    figure; plot(time(2:end),series)
    hold on
    plot(time,thres*ones(length(time),1),'k--')
    if ~isempty(clusters_start)
        plot((clusters_start/fs+presec),thres,'g.','markersize',10)
        plot((clusters_end/fs+presec),thres,'r.','markersize',10)
    end
    plot([3600,3600],[min(series),max(series)])
end