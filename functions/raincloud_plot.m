function [h, u] = raincloud_plot(X, varargin)
% RAINCLOUD_PLOT - Violin + box + raw data plot (Micah Allen style)
% Skips if X is empty or NaN. Shows n= count. Clean layout.
% -------------------------------------------------------------------------

% Remove NaNs and skip if empty
X = X(~isnan(X));
if isempty(X)
    h = {};
    u = NaN;
    return;
end

% -------------------- INPUT PARSING --------------------
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p, 'X', @isnumeric);
addOptional(p, 'color', [0.5 0.5 0.5], @isnumeric)
addOptional(p, 'band_width', [])
addOptional(p, 'density_type', 'ks', @ischar)
addOptional(p, 'box_on', 0, @isnumeric)
addOptional(p, 'box_dodge', 0, @isnumeric)
addOptional(p, 'box_dodge_amount', 0.5, @isnumeric)
addOptional(p, 'alpha', 1, validScalarPosNum)
addOptional(p, 'dot_dodge_amount', 0.6, @isnumeric)
addOptional(p, 'box_col_match', 1, @isnumeric)
addOptional(p, 'line_width', 2, validScalarPosNum)
addOptional(p, 'lwr_bnd', 1, @isnumeric)
addOptional(p, 'bxcl', [0 0 0], @isnumeric)
addOptional(p, 'bxfacecl', [1 1 1], @isnumeric)
addOptional(p, 'cloud_edge_col', 'none', @(x) (ischar(x) && isstringlike(x)) || (isnumeric(x) && numel(x) == 3))
addOptional(p, 'group_label', '', @ischar)
addOptional(p, 'x_limits', [], @(x) isnumeric(x) && numel(x) == 2);

parse(p, X, varargin{:});

% -------------------- ASSIGN PARSED INPUT --------------------
color               = p.Results.color;
density_type        = p.Results.density_type;
box_on              = p.Results.box_on;
box_dodge           = p.Results.box_dodge;
box_dodge_amount    = p.Results.box_dodge_amount;
alpha               = p.Results.alpha;
dot_dodge_amount    = p.Results.dot_dodge_amount;
box_col_match       = p.Results.box_col_match;
line_width          = p.Results.line_width;
lwr_bnd             = p.Results.lwr_bnd;
bxcl                = p.Results.bxcl;
bxfacecl            = p.Results.bxfacecl;
cloud_edge_col      = p.Results.cloud_edge_col;
band_width          = p.Results.band_width;
group_label         = p.Results.group_label;
x_limits            = p.Results.x_limits;

% If empty, default to color
if isempty(cloud_edge_col)
    cloud_edge_col = color;
end

% -------------------- KERNEL DENSITY ESTIMATION --------------------
switch density_type
    case 'ks'
        [f, Xi, u] = ksdensity(X, 'bandwidth', band_width, 'NumPoints', 200);
    case 'rash'
        assert(exist('rst_RASH', 'file') == 2, ...
            'RASH requires Robust Stats Toolbox');
        [Xi, f] = rst_RASH(X);
        u = NaN;
end

% -------------------- LAYOUT --------------------
hold on;
max_f = max(f);
wdth = max_f * 0.25;

cloud_offset = 0;
box_offset = -max_f * 0.3;
drop_offset = -max_f * 0.8;

% -------------------- CLOUD PLOT FIRST --------------------
h{1} = area(Xi, f + cloud_offset, 'FaceColor', color, ...
    'EdgeColor', cloud_edge_col, 'LineWidth', line_width, 'FaceAlpha', alpha);

% -------------------- BOX PLOT --------------------
quartiles = quantile(X, [0.25 0.75 0.5]);
iqr = quartiles(2) - quartiles(1);
Xs = sort(X);

whiskers = [quartiles(1), quartiles(2)];
low_whisker = Xs(Xs > (quartiles(1) - 1.5 * iqr));
high_whisker = Xs(Xs < (quartiles(2) + 1.5 * iqr));
if ~isempty(low_whisker), whiskers(1) = min(low_whisker); end
if ~isempty(high_whisker), whiskers(2) = max(high_whisker); end

if box_on
    if box_col_match, bxcl = color; end
    box_height = wdth * 1.0;

    % Box
    h{2} = rectangle('Position', [quartiles(1), box_offset - box_height/2, ...
        quartiles(2)-quartiles(1), box_height], ...
        'EdgeColor', bxcl, 'LineWidth', line_width, 'FaceColor', bxfacecl);

    % Median
    h{3} = line([quartiles(3) quartiles(3)], ...
        [box_offset - box_height/2, box_offset + box_height/2], ...
        'Color', bxcl, 'LineWidth', line_width);

    % Whiskers
    h{4} = line([whiskers(1) quartiles(1)], [box_offset box_offset], ...
        'Color', bxcl, 'LineWidth', line_width);
    h{5} = line([quartiles(2) whiskers(2)], [box_offset box_offset], ...
        'Color', bxcl, 'LineWidth', line_width);
end

% -------------------- SCATTER (DROPLETS) --------------------
jit = (rand(size(X)) - 0.5) * wdth;
drops_pos = drop_offset + jit;
h{6} = scatter(X, drops_pos, 10, ...
    'MarkerFaceColor', color, 'MarkerEdgeColor', 'none');

% -------------------- n= TEXT --------------------
% Clamp inside plot
if isempty(x_limits)
    xlims = [min(X) max(X)];
else
    xlims = x_limits;
end

x_text = xlims(2) - 0.02 * range(xlims);
y_text = cloud_offset + max(f) * 0.5;
n_text = sprintf('n = %d', numel(X));
h{7} = text(x_text, y_text, n_text, ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 11, 'FontWeight', 'bold');

% -------------------- GROUP LABEL --------------------
if ~isempty(group_label)
    title(group_label, 'Interpreter', 'none');
end

% -------------------- LIMITS --------------------
xlim(xlims);
ylim([-max_f * 1.3, max(f) * 1.2]);
end