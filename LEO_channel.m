%% Rayleigh Channel example using NTN-TDL-A channel model
max_speed = 500; % unit km/h
delays = [0 1.0811 2.8416]*30*10^(-9); % multiply delay spread 30ns and delay unit is ns.
pdp = [0 -4.675 -6.482]; % power delay profile (dB)
taps = length(delays); % number of delay taps
fc = 2.55*10^9; % carrier frequency=2.55Ghz
fs = 61.44*10^6; % sample rate=61.44Mhz
c = physconst('LightSpeed');
max_Doppler = fc*(max_speed*1000/3600)/c;

rayleighchan = comm.RayleighChannel(...  
    'SampleRate', fs,...
    'PathDelays', delays,...
    'AveragePathGains', pdp,...
    'MaximumDopplerShift', max_Doppler,...
    'DopplerSpectrum', doppler('Jakes')...
    );
Rx_waveform = rayleighchan(Tx_waveform);
Rx_waveform = awgn(Rx_waveform, SNR, 'measured');


%% Rician Channel example using NTN-TDL-D channel model
max_speed = 500; % unit km/h
delays = [0 0.5596 7.3340]*30*10^(-9);  % multiply delay spread 30ns and delay unit is ns.
pdp = [-0.284 -9.887 -16.771]; % power delay profile (dB)
taps = length(delays); % number of delay taps
fc = 2.55*10^9; % carrier frequency=2.55Ghz
fs = 61.44*10^6; % sample rate=61.44Mhz
c = physconst('LightSpeed');
max_Doppler = fc*(max_speed*1000/3600)/c;
KFactor = 11.707; % unit dB

ricianchan = comm.RicianChannel(...  
    'SampleRate', fs,...
    'PathDelays', delays,...
    'AveragePathGains', pdp,...
    'KFactor', db2pow(KFactor),...    % K-Factor要代linear power ratio
    'DirectPathDopplerShift', max_Doppler,...
    'MaximumDopplerShift', max_Doppler,...
    'DopplerSpectrum', doppler('Jakes')...
    );
Rx_waveform = ricianchan(Tx_waveform);
Rx_waveform = awgn(Rx_waveform, SNR, 'measured');