% Simulate starting with FFT -> then convert to time-domain signal. 

%% 

clear
addpath(genpath('lib'))

fs = 100; 

dur = 2.4 * 3; 

N = round(fs * dur); 

mod(N, 2)

hN = N / 2 + 1; 


frex = 1/2.4 * [1 : 12]; 
frex_idx = round(frex / fs * N + 1); 

magnitudes = zeros(1, N); 
phases = zeros(1, N); 

magnitudes(frex_idx) = rand(1, length(frex)); 
phases(frex_idx) = 2 * pi * rand(1, length(frex)); 

% negative frequencies
magnitudes(hN+1 : end) = flip(magnitudes(2 : hN - 1)); 
phases(hN+1 : end) = -flip(phases(2 : hN - 1)); 

X = magnitudes .* exp(1j .* phases); 

freq = [0 : N-1] / N * fs; 


x = ifft(X); 

isreal(x)



f = figure('color', 'white', 'pos', [297 363 1215 419]); 
pnl = panel(f); 
pnl.pack('h', [70, 30]); 

pnl(1).pack('v', 2); 
pnl(2).pack({[0, 0, 1, 1]}); 
pnl(2).pack({[0, 0.3, 1, 0.3]}); 

pnl(1, 1).select(); 
stem(freq, magnitudes)
pnl(1, 1).ylabel('magnitude'); 

pnl(1, 2).select(); 
stem(freq, phases)
pnl(1, 2).ylabel('phase (radians)'); 

pnl(2, 2).select(); 
t = [0 : N-1] / fs; 
plot(t, x, 'linew', 2); 


pnl(1).xlabel('frequency (Hz)'); 
pnl(2, 2).xlabel('time (s)'); 

pnl.fontsize = 12; 


%% from reply to Rajendran PNAS 2018

clear

fs = 500; 

pattern = [1 1 1 0 1 1 1 0 1 1 0 0]; 

IOI = 0.2; 
rampon = 0.01; 
rampoff = 0.05; 

ncycles = 24; 

t = [0 : round(ncycles*length(pattern)*IOI*fs)-1]/fs; 

% create the stimulus 
env = zeros(1, length(t)); 
env_event = ones(1,round(IOI*fs)); 
env_event(1:round(rampon*fs)) = linspace(0,1,round(rampon*fs)); 
env_event(end-round(rampoff*fs)+1:end) = linspace(1,0,round(rampoff*fs)); 
for i=0:ncycles*length(pattern)-1
    if pattern(mod(i,length(pattern))+1)        
        env(round(i*IOI*fs)+1:round(i*IOI*fs+length(env_event))) = env_event; 
    end
end

N = length(env); 
hN = floor(N/2)+1; 
X = fft(env);

idx5Hz = round(5/fs*N)+1; 
n_beats = length(t)/fs/(IOI*4); 
beat_times = IOI*4 * [1:n_beats]; 
pattern_times = IOI*12 * [1:ncycles]; 

% original (0 to 5 Hz)
Xshifted = fftshift(X); 
Xshifted([1:hN-idx5Hz-1, hN+idx5Hz:end]) = 0; 
env_original = ifft(ifftshift(Xshifted)); 

% zero phase
mX = abs(Xshifted); 
X_zeroPhase = mX * exp(0i); 
env_zeroPhase = ifft(ifftshift(X_zeroPhase)); 

% random phase
phases_plus = (rand(1,N/2-1)*2*pi-pi); 
phases_minus = fliplr(-phases_plus); 
phases = [0, phases_minus, 0, phases_plus]; 
X_randPhase = mX .* exp(phases*1i); 
env_randPhase = ifft(ifftshift(X_randPhase)); 


% -------------------- PLOT --------------------
figure('color','white')
subplot 311
plot(t,env_original,'k','LineWidth',1.7); 
hold on
plot(beat_times,max(env_original)*1.1,'ro','MarkerFaceColor','r'); 
xlim([0,2.4*3])
box off

subplot 312
plot(t,env_zeroPhase,'LineWidth',1.7,'color',[0.4902,0.1804,0.5608]); 
hold on
plot(beat_times,max(env_zeroPhase)*1.1,'ro','MarkerFaceColor','r'); 
xlim([0,2.4*3])
box off

subplot 313
plot(t,env_randPhase,'LineWidth',1.7,'color',[0.4706,0.6706,0.1882]); 
hold on
plot(beat_times,max(env_randPhase)*1.1,'ro','MarkerFaceColor','r'); 
xlim([0,2.4*3])
box off
% -------------------------------------------------

