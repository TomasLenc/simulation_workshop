% Illustrate how FFT is computed just using cosine waves. 
clear
addpath(genpath('lib'))

N = 32; 
n = [0 : N-1]; 

col = [50, 168, 82]/255; 
linew = 3; 

%% generating sine waves

f = figure('color', 'white', 'pos', [637 556 475 400]); 
pnl = panel(f); 
pnl.pack('v', 2); 


ax = pnl(1).select(); 

omg = (2 * pi * 3) / N; 
x = cos(omg * n); 
stem(n, x, 'color', col, 'linew', linew, 'marker', 'o'); 
ax.YLim = [-1, 1]; 
ax.YTick = [-1, 1]; 
ax.XLim = [0, N-1]; 

ax = pnl(2).select(); 

A = 3; 
phi = pi / 3; 
omg = (2 * pi * 3) / N; 
x = A * cos(omg * n - phi); 
stem(n, x, 'color', col, 'linew', linew, 'marker', 'o'); 
ax.YLim = [-3, 3]; 
ax.YTick = [-3, 3]; 
ax.XLim = [0, N-1]; 

pnl.fontsize = 12; 
pnl.xlabel('sample');
pnl.de.margin = [10, 15, 5, 5]; 



%% sinewave with any A and phi is just a weighted cos + sin

f = figure('color', 'white', 'pos', [632 564 544 380]); 
pnl = panel(f); 
pnl.pack('v', 2); 
pnl(1).pack('h', [70, 30]); 
pnl(2).pack('h', [70, 30]); 

A = 3; 
phi = pi / 3; 
omg = (2 * pi * 3) / N; 

x = A * cos(omg * n - phi); 

ax = pnl(1, 1).select(); 
plot(n, x, 'color', 'black', 'linew', linew); 
ax.YLim = [-3, 3]; 
ax.YTick = [-3, 3]; 
ax.XLim = [0, N-1]; 

ax = polaraxes; 
polarplot([phi, phi], [0, A], 'linew', linew); 
ax.RTick = [3]; 
ax.FontSize = 12; 
pnl(1, 2).select(ax); 

c = A * cos(phi); 
cos_wave = cos(omg * n); 

s = A * sin(phi); 
sin_wave = sin(omg * n); 

x_sum = c * cos_wave + s * sin_wave; 

ax = pnl(2, 1).select(); 
h = plot(n, c * cos_wave, 'linew', linew); 
h.Color(4) = 0.5; 
hold on
h = plot(n, s * sin_wave, 'linew', linew); 
h.Color(4) = 0.5; 
plot(n, x_sum, 'color', 'black', 'linew', linew); 
ax.YLim = [-3, 3]; 
ax.YTick = [-3, 3]; 
ax.XLim = [0, N-1]; 

ax = pnl(2, 2).select(); 
plot([0, c], [0, s], 'linew', linew); 
ax.XLim = [-4, 4]; 
ax.YLim = [-4, 4]; 
ax.XTick = [-3, 3]; 
ax.YTick = [-3, 3]; 
xlabel('c'); 
ylabel('s'); 

ax.XAxisLocation = 'origin'; 
ax.YAxisLocation = 'origin'; 

pnl.fontsize = 12; 
pnl(2, 1).xlabel('sample');
pnl.de.margin = [10, 15, 5, 5]; 
pnl.marginright = 15; 


%% exponential numbers are just a convenient way to do arithmetics with cos and sin

% Euler's formula
%     exp(i*phi) = cos(phi) + i*sin(phi)

A = 3; 
phi = pi / 3; 

c = A * cos(phi); 
cos_wave = cos(omg * n); 

s = A * sin(phi); 
sin_wave = sin(omg * n); 


x = c + 1j*s; 

y = A * exp(1j * phi); 


% easy to show how multiplying 2 sinewaves works
% 
%     A * exp(i*alpha)   *   B * exp(i*beta) =  A * B * exp(i*[alpha+beta])
% 
% -> exp representation super handy if you want to change phase


%% DFT matrix

linew = 2; 

omg = (2 * pi * 0) / N; 
col0_c = cos(omg * n); 
col0_s = sin(omg * n); 

omg = (2 * pi * 1) / N; 
col1_c = cos(omg * n); 
col1_s = sin(omg * n); 

omg = (2 * pi * 2) / N; 
col2_c = cos(omg * n); 
col2_s = sin(omg * n); 

omg = (2 * pi * 3) / N; 
col3_c = cos(omg * n); 
col3_s = sin(omg * n); 


f = figure('color', 'white', 'pos', [98 139 1747 578]); 
pnl = panel(f); 
pnl.pack('h', 16+1); 

cols = brewermap(16+1, '-Set3'); 

for i_col=0:16
    
    pnl(i_col+1).pack('h', 2); 

    omg = (2 * pi * i_col) / N; 
    col_c = cos(omg * n); 
    col_s = sin(omg * n); 

    ax = pnl(i_col+1, 1).select(); 
    stem(n, col_c, 'linew', linew, 'color', cols(i_col+1, :)); 
    ax.XLim = [0, N-1];
    ax.XAxis.Visible = 'off'; 
    ax.YAxis.Visible = 'off'; 
    view(ax, 90, 90); 
    
    ax = pnl(i_col+1, 2).select(); 
    stem(n, col_s, 'linew', linew, 'color', cols(i_col+1, :)); 
    ax.XLim = [0, N-1];
    ax.YLim = [-1, 1]; 
    ax.XAxis.Visible = 'off'; 
    ax.YAxis.Visible = 'off'; 
    view(ax, 90, 90); 
    
    pnl(i_col+1).title(sprintf('%d', i_col))
    
end

pnl.de.margin = 0; 
pnl.de.marginright = 3; 

pnl.margintop = 15; 
pnl.fontsize = 12; 



%% DFT sinewaves in polar coordinates 
% + negative frequencies? 

linew = 2; 

% try with 3, -3, 32-3
n_cycles = 32 -3; 


omg = (2 * pi * n_cycles) / N; 


f = figure('color', 'white', 'pos', [98 139 800 300]); 
pnl = panel(f); 

pnl.pack('h', [70, 30]); 
pnl(1).pack('v', 2); 


ax_cos = pnl(1, 1).select(); 
ax_cos.XLim = [0, N-1]; 
ax_cos.YLim = [-1, 1]; 
stem(0, cos(omg * 0), 'o', ...
    'color', 'r', 'markerfacecolor', 'r', 'markersize', 10); 
hold(ax_cos, 'on'); 


ax_sin = pnl(1, 2).select(); 
ax_sin.XLim = [0, N-1]; 
ax_sin.YLim = [-1, 1]; 
plot(0, sin(omg * 0), 'o', ...
    'color', 'b', 'markerfacecolor', 'b', 'markersize', 10); 
hold(ax_sin, 'on'); 



ax_polar = polaraxes; 
hold on
polarplot([pi, 0], [1, 1], 'linew', 3, 'color', 'k')
polarplot([-pi/2, pi/2], [1, 1], 'linew', 3, 'color', 'k')

ax_polar.RTick = []; 
ax_polar.RLim = [0, 1]; 
ax_polar.ThetaAxisUnits = 'radians'; 
ax_polar.ThetaTick = 2*pi * [0:N-1]/N; 
ax_polar.ThetaTickLabel = cellfun(@num2str, num2cell([0:N-1]), 'uni', 0); 

pnl(2).select(ax_polar); 
h_polar = polarplot(omg * 0, 1, 'o',...
    'color', 'k', 'markerfacecolor', 'k', 'markersize', 10); 



% plot thre rest of the points one by one
for n=1:N-1
    
    h_polar.ThetaData = omg * n; 
    
    stem(ax_cos, n, cos(omg * n), 'o', ...
        'color', 'r', 'markerfacecolor', 'r', 'markersize', 10); 
    stem(ax_sin, n, sin(omg * n), 'o', ...
        'color', 'b', 'markerfacecolor', 'b', 'markersize', 10); 
    
    pause(0.1); 
end

% plot everything together
delete(h_polar); 

polarplot(ax_polar, omg * [0:N-1], 0.95, 'o', ...
    'color', 'k', 'markerfacecolor', 'k', 'markersize', 10); 

pnl.margin = [10, 10, 10, 10]; 

pnl.title(sprintf('%d cycles per vector', n_cycles)); 
pnl.fontsize = 12; 








