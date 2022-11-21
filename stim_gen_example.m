% Here we use convolution to generate signals composed of events with a
% specific shape, arranged in time to create a particular temporal pattern. 

% Check out documentation for the rnb_tools functions for more :) 

clear 
fprintf('\n');
addpath(genpath('lib'))
addpath(genpath('/home/tomo/Documents/MATLAB/rnb_tools/src')); 

%% simulate

fs = 200; 

pat = [1 0 1 1 1 1 0 1 1 1 0 0] % [1 1 1 0 1 1 1 0 1 1 0 0], % [1 0 1 1 1 1 0 1 1 1 0 0]
    
grid_ioi = 0.2; 

ir_type = 'erp2'; 

emph_levels = linspace(0, 1, 6); 

n_cond = length(emph_levels); 

% colors
cmap_name = '-YlGn'; 
colors = num2cell(brewermap(n_cond + n_cond, cmap_name), 2); 
colors = colors(1:n_cond, :); 

fontsize = 14; 

%% 

if strcmp(ir_type, 'square')
    ir = get_square_kernel(fs, ...
        'duration', 0.100, ...
        'rampon', 0, ...
        'rampoff', 0 ...
        ); 
elseif strcmp(ir_type, 'erp')
    ir = get_erp_kernel(fs,...
        'amplitudes', 1,...
        't0s', 0, ...
        'taus', 0.050, ...
        'f0s', 7, ...
        'duration', 0.2 ...
        ); 
elseif strcmp(ir_type, 'erp2')
    ir = get_erp_kernel(fs,...
        'amplitudes', [0.4, 0.75],...
        't0s', [0, 0], ...
        'taus', [0.2, 0.050], ...
        'f0s', [1, 7], ...
        'duration', 0.5 ...
        ); 
end

f = figure('color','white', 'position',[575 1 582 500]); 
pnl = panel(f); 
pnl.pack('v', n_cond + 1); 
for i=1:n_cond+1
    pnl(i).pack('h', [45, 10, 45]); 
end

ax = pnl(1, 1).select(); 
t_pattern = [0 : length(pat)-1] * 0.2; 
scatter(t_pattern(logical(pat)), zeros(1, sum(pat)), 100, 'kx', 'LineWidth', 2)
hold on
scatter(t_pattern(logical(~pat)), zeros(1, sum(~pat)), 50, 'k.', 'LineWidth', 2)
ax.XLim = [0, 2.4]; 
ax.Visible = 'off'; 

ax = pnl(4, 2).select(); 
t_ir = [0 : length(ir)-1] / fs; 
plot(t_ir, ir, 'k', 'linew', 3); 
ax.YAxis.Visible = 'off'; 
ax.XTick = [0, max(ax.XTick)]; 
ax.TickDir = 'out'; 

for i_cond=1:n_cond
    
    emph = emph_levels(i_cond); 

    fprintf('calculating %d/%d\n', i_cond, n_cond)
    
    % make whole signal 
    [x, t, x_dirac] = get_s(...
                        pat, ...
                        0.2, ...
                        fs, ...
                        'n_cycles', 16, ...
                        'ir', ir, ...
                        'emph_magn', emph, ...
                        'emph_period', 4 ...
                        );
                
    ax = pnl(i_cond + 1, 1).select(); 
    stem(t, x_dirac, 'linew', 3, 'marker', 'none', 'color', colors{i_cond}); 
    ax.XLim = [0, 2.4-1/fs]; 
    ax.Visible = 'off'; 
    
    ax = pnl(i_cond + 1, 3).select(); 
    plot(t, x, 'linew', 3, 'color', colors{i_cond}); 
    ax.XLim = [0, 2.4-1/fs]; 
    ax.Visible = 'off'; 
    
end

pnl.de.margin = 5; 

pnl(1).marginbottom = 0; 

pnl.margin = [5, 5, 5, 5]; 

pnl.fontsize = 14; 

