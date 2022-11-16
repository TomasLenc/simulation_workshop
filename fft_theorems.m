
% https://ccrma.stanford.edu/~jos/

% http://www.cns.nyu.edu/~eero/math-tools19/


%% shift - multiplication

clear
addpath(genpath('lib'))

fs = 500; 
dur = 10; 
N = round(fs * dur); 
hN = floor(N / 2) + 1; 
t = [0 : N-1] / fs; 
freq = [0 : hN-1] / N * fs; 

col1 = [50, 168, 82]/255; 
col2 = [168, 144, 50]/255; 
col = [0, 0, 0]/255; 

x1 = sin(2 * pi * t * 100); 
x2 = 1 * sin(2 * pi * t * 3) + ...
     0.5 * sin(2 * pi * t * 11) + ...
     2 * sin(2 * pi * t * 0) ; 
x = x1 .* x2; 

mX1 = abs(fft(x1)); 
mX2 = abs(fft(x2)); 
mX = abs(fft(x)); 


f = figure('color', 'white', 'pos',[637 370 969 387]); 
pnl = panel(f); 

pnl.pack('h', 2); 
pnl(1).pack('v', 3); 
pnl(2).pack('v', 3); 

linew = 2; 

ax = pnl(1, 1).select(); 
plot(t, x1, 'color', col1, 'linew', linew); 
ax.XLim = [0, 1]; 
ax.XTick = []; 

ax = pnl(1, 2).select(); 
plot(t, x2, 'color', col2, 'linew', 3); 
ax.XLim = [0, 1]; 
ax.XTick = []; 

ax = pnl(1, 3).select(); 
plot(t, x, 'color', col, 'linew', linew); 
ax.XLim = [0, 1]; 


ax = pnl(2, 1).select(); 
stem(freq, mX1(1:hN), 'color', col1, 'linew', linew, 'marker', 'none'); 
ax.XLim = [0, fs/2]; 
ax.YTick = []; 
ax.XTick = []; 

ax = pnl(2, 2).select(); 
stem(freq, mX2(1:hN), 'color', col2, 'linew', linew, 'marker', 'none'); 
ax.XLim = [0, fs/2]; 
ax.YTick = []; 
ax.XTick = []; 

ax = pnl(2, 3).select(); 
stem(freq, mX(1:hN), 'color', col, 'linew', linew, 'marker', 'none'); 
ax.XLim = [0, fs/2]; 
ax.YTick = []; 

pnl.fontsize = 12; 
pnl(1).xlabel('time (s)');
pnl(2).xlabel('freq (Hz)');

pnl.de.margin = [10, 5, 5, 5]; 


%% zero-padding - approaching continuous fourier in the limit

% This about what we're doing: 
%     1) window infinite sine wave with a rectangular window
%     2) sample the continuous FT

clear
addpath(genpath('lib'))

fs = 500; 
dur = 0.5; 
N = round(fs * dur); 
hN = floor(N / 2) + 1; 
t = [0 : N-1] / fs; 
freq = [0 : hN-1] / N * fs; 

N_pad = N + N*5; 
hN_pad = floor(N_pad/2) + 1; 
t_pad = [0 : N_pad-1] / N_pad * fs; 
freq_pad = [0 : hN_pad-1] / N_pad * fs; 

col = [50, 168, 82]/255; 
col_pad = [168, 144, 50]/255; 

x = sin(2 * pi * t * 10); 
mX = abs(fft(x)); 

x_pad = [x, zeros(1, N_pad -N)]; 
mX_pad = abs(fft(x_pad)); 


f = figure('color', 'white', 'pos',[637 370 969 387]); 
pnl = panel(f); 

pnl.pack('h', 2); 
pnl(1).pack('v', 2); 
pnl(2).pack('v', 2); 

linew = 2; 

ax = pnl(1, 1).select(); 
plot(t, x, 'color', col, 'linew', linew); 

ax = pnl(1, 2).select(); 
plot(t_pad, x_pad, 'color', col_pad, 'linew', 3); 

ax = pnl(2, 1).select(); 
stem(freq, mX(1:hN), 'color', col, 'linew', linew, 'marker', 'none'); 
ax.XLim = [0, fs/2]; 
ax.YTick = []; 
ax.XTick = []; 
ax.XLim = [0, 40]; 

ax = pnl(2, 2).select(); 
plot(freq_pad, mX_pad(1:hN_pad), 'color', col_pad, 'linew', linew); 
ax.XLim = [0, fs/2]; 
ax.YTick = []; 
ax.XTick = []; 
ax.XLim = [0, 40]; 

pnl.fontsize = 12; 
pnl(1).xlabel('time (s)');
pnl(2).xlabel('freq (Hz)');

pnl.de.margin = [10, 5, 5, 5]; 



%% stretch - repeat

clear
addpath(genpath('lib'))
addpath(genpath('/home/tomo/Documents/MATLAB/rnb_tools'))

fs = 100; 

col1 = [50, 168, 82]/255; 
col2 = [168, 144, 50]/255; 
col = [0, 0, 0]/255; 

nreps = [1 : 4]; 

x_unit = get_colored_noise(round(1 * fs), fs, -2); 

f = figure('color', 'white', 'pos',[637 370 969 387]); 
pnl = panel(f); 

pnl.pack('h', 2); 
pnl(1).pack('v', length(nreps)); 
pnl(2).pack('v', length(nreps)); 

linew = 2; 

max_dur = (max(nreps) * length(x_unit)) / fs; 

for i=1:length(nreps)
    
    nrep = nreps(i); 
    
    x = repmat(x_unit, 1, nrep); 
    mX = abs(fft(x)); 

    N = length(x); 
    hN = floor(N / 2) + 1; 
    t = [0 : N-1] / fs; 
    freq = [0 : hN-1] / N * fs; 

    ax = pnl(1, i).select(); 
    plot(t, x, 'color', col1, 'linew', linew); 
    ax.XLim = [0, max_dur]; 
    ax.XTick = [0, max_dur]; 
    ax.YTick = []; 

    ax = pnl(2, i).select(); 
    stem(freq, mX(1:hN), 'color', col1, 'linew', linew, 'marker', 'o'); 
    ax.XLim = [0, fs/10]; 
    ax.XTick = [0, fs/10]; 
    ax.YTick = []; 

end

pnl.fontsize = 12; 
pnl(1).xlabel('time (s)');
pnl(2).xlabel('freq (Hz)');

pnl.de.margin = [10, 8, 5, 5]; 



%% convolution theorem

clear
addpath(genpath('lib'))
addpath(genpath('/home/tomo/Documents/MATLAB/rnb_tools'))

fs = 100; 
dur = 0.2 * 12 * 4; 
N = round(fs * dur); 
hN = floor(N / 2) + 1; 
t = [0 : N-1] / fs; 
freq = [0 : hN-1] / N * fs; 

x = get_s([1,1,1,1,0,1,1,1,0,0,1,0], 0.2, fs, 'n_cycles', 4); 
mX = abs(fft(x)); 

ir = get_square_kernel(fs, 'duration', 0.100); 
t_ir = [0 : length(ir)-1] / fs; 

ir_pad = [ir, zeros(1, N-length(ir))]; 
mX_ir = abs(fft(ir_pad)); 

x_conv_ir = conv(x, ir, 'same'); 
mX_s_conv_ir = abs(fft(x_conv_ir)); 


col1 = [50, 168, 82]/255; 
col2 = [168, 144, 50]/255; 
col = [0, 0, 0]/255; 


f = figure('color', 'white', 'pos',[637 370 969 387]); 
pnl = panel(f); 

pnl.pack('h', 2); 
pnl(1).pack('v', 3); 
pnl(2).pack('v', 3); 

linew = 2; 

ax = pnl(1, 1).select(); 
plot(t, x, 'color', col1, 'linew', linew); 

pnl(1, 2).pack({[0, 0, 1, 1]}); 
pnl(1, 2).pack({[0.7, 0.5, 0.3, 0.6]}); 

ax = pnl(1, 2, 1).select(); 
plot(t, ir_pad, 'color', col2, 'linew', 3); 

ax = pnl(1, 2, 2).select(); 
plot(t_ir, ir, 'color', col2, 'linew', 3); 
ax.YTick = []; 

ax = pnl(1, 3).select(); 
plot(t, x_conv_ir, 'color', col, 'linew', linew); 


ax = pnl(2, 1).select(); 
stem(freq, mX(1:hN), 'color', col1, 'linew', linew, 'marker', 'none'); 
ax.XLim = [0, fs/2]; 
ax.YTick = []; 
ax.XTick = []; 

ax = pnl(2, 2).select(); 
stem(freq, mX_ir(1:hN), 'color', col2, 'linew', linew, 'marker', 'none'); 
ax.XLim = [0, fs/2]; 
ax.YTick = []; 
ax.XTick = []; 

ax = pnl(2, 3).select(); 
stem(freq, mX_s_conv_ir(1:hN), 'color', col, 'linew', linew, 'marker', 'none'); 
ax.XLim = [0, fs/2]; 
ax.YTick = []; 

pnl.fontsize = 12; 
pnl(1).xlabel('time (s)');
pnl(2).xlabel('freq (Hz)');

pnl.de.margin = [10, 5, 5, 5]; 











