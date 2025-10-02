function [found1, found2] = check_signals_dual(x, time, events, t_limits)
% CHECK_SIGNALS_DUAL - Visual signal quality assessment for mean and std signals
%
% Inputs:
%   x         - 3D array (time x trials x 2), x(:,:,1)=mean, x(:,:,2)=std
%   time      - time vector
%   events    - event stop times vector (e.g. one stop time per trial)
%   t_limits  - x-axis limits, e.g. [-12, 12]
%
% Outputs:
%   found1    - logical vector (1=keep, 0=discard) for mean
%   found2    - logical vector for std

dbg_fig = 100;
n_trials = size(x, 2);
found1 = true(1, n_trials);
found2 = true(1, n_trials);

% Create figure and plot all signals
figure(dbg_fig);
clf;
r = ceil(sqrt(n_trials));
c = ceil(n_trials / r);
set(gcf, 'Position', [100, 100, 1200, 800]);
sgtitle('Click plots to mark bad signals (yellow background). Press any key when done.');

for k = 1:n_trials
    subplot(r, c, k);
    
    yyaxis left
    plot(time, x(:, k, 1), 'b'); % mean
    ylabel('Mean');

    yyaxis right
    plot(time, x(:, k, 2), 'r'); % std
    ylabel('STD');

    hold on;

    g = gca;

    % Draw time zero marker
    plot([0, 0], [g.YLim(1), g.YLim(2)], 'k--', 'linewidth', 2);

    % Draw t_stop for this trial only
    plot([events(k), events(k)], [g.YLim(1), g.YLim(2)], 'k--', 'linewidth', 2);

    xlim(t_limits);

    set(gca, 'ButtonDownFcn', {@toggle_axes, k});
    ch = get(gca, 'Children');
    for i = 1:numel(ch)
        set(ch(i), 'ButtonDownFcn', {@toggle_axes, k});
    end
end

set(findall(gcf, 'type', 'axes'), 'FontSize', 14, 'FontName', 'Times');

% Wait for user interaction
disp('Click plots to mark bad signals. Press any key when done.');
pause;

% Collect signal quality flags
ax_all = findall(gcf, 'type', 'axes');
for k = 1:n_trials
    ax = ax_all(end - k + 1); % subplot order
    is_bad = isequal(get(ax, 'Color'), [1, 1, 0]);
    found1(k) = ~is_bad;
    found2(k) = ~is_bad; % both use same logic
end

close(gcf);

end

function toggle_axes(src, ~, idx)
    if strcmp(get(src, 'Type'), 'axes')
        ax = src;
    else
        ax = get(src, 'Parent');
    end
    current_color = get(ax, 'Color');
    if isequal(current_color, [1, 1, 0])
        set(ax, 'Color', 'w');
    else
        set(ax, 'Color', [1, 1, 0]);
    end
end
