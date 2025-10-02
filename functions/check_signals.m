function found = check_signals(x, time, events)
% - function found = check_signals(tld) -- visualises 
% time series for user to decide which signals should be kept and discarded
%
%
% Input
% - x              - Time-locked heart rate, saturation, and respiratory
%                    rate data as double variable
% - time           - Corresponding time vector as double variable
%
% Output
% - found          - Signal quality (0: discarded; 1: include) as double
%                    variable
%
%
% See also MOVING_AVERAGE, TIME_LOCK_DATA, PLOT_VITAL_SIGNS, FIND_RESPONDERS
% ________________________________________________________________________
%
% This file is released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html
%
%                                           (c) Coen Zandvoort, 2024
% ________________________________________________________________________

if nargin == 1 && ishandle(x) % function is callback after button-press
    switch lower(get(x, 'Type'))
        case 'axes', ax = x;
        case 'line', ax = get(x, 'Parent');
    end
    set(ax, 'color', xor(get(ax, 'Color'), [0, 0, 1])); % toggle the bg-color
    found = updateSignalPlot(get(ax, 'Parent'));
    return
end


dbg_fig = 100; % figure number
found = true(1, numel(x));


% graphical debugging or interactive mode
plot_signals(dbg_fig(1), time, x, events, found);
set(findobj(dbg_fig(1), 'type', 'axes'), 'buttondownfcn', [mfilename '(gcbo);']);
set(findobj(dbg_fig(1), 'type', 'line'), 'buttondownfcn', [mfilename '(gcbo);']);

fprintf('use the mouse to select signals and press any key to finish selection ');
pause; % wait for key-press
fprintf('\n')
found = updateSignalPlot(dbg_fig(1));

found = found(:);


% plotting signals
function plot_signals(fig, t, x, ev, found)
Mm = size(x, 2);
r = ceil(sqrt(Mm));
u.ax = [];
figure(fig); 
po = get(gcf, 'position');
set(gcf, 'position', [po(1 : 2), 900, 600])
clf;
for k = 1 : Mm
    subplot(r, ceil(Mm / r), k);
    plot(t, x(:, k));
    hold on

    g = gca;
    for e = 1 : size(ev, 2)
        plot([0, 0], [g.YLim(1), g.YLim(2)], 'k--', 'linewidth', 2)
        plot([ev(e, 1), ev(e, 1)], [g.YLim(1), g.YLim(2)], 'k--', 'linewidth', 2)
    end

    u.ax(k) = gca;
    set(u.ax(~found(k)), 'color', [1, 1, 0]);
    set(u.ax, 'xlim', t([1, end]));
end

set(fig, 'name', 'signals', 'UserData', u);
set(findall(0, 'type', 'axes'), 'Fontsize', 16, 'FontName', 'Times');


% update signals (callback function when clicking on signals)
function found = updateSignalPlot(figModes)

ax = getfield(get(figModes, 'UserData'), 'ax');
found = true(size(ax));
if ~isempty(ax)
    
    h = get(ax, 'color');
    for i = 1 : numel(h)
        found(i) = ~(sum(h{i} == [1, 1, 0]) == 3);
    end

end

% _ EOF____________________________________________________________________