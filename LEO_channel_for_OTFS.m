%% Rayleigh Channel example using NTN-TDL-A channel model
max_speed = 500; % unit km/h
delays = [0 1.0811 2.8416]*30*10^(-9); % multiply delay spread 30ns and delay unit is ns.
pdp = [0 -4.675 -6.482]; % power delay profile (dB)
taps = length(delays); % number of delay taps
fc = 2.55*10^9; % carrier frequency=2.55Ghz
fs = 61.44*10^6; % sample rate=61.44Mhz
c = physconst('LightSpeed');
max_Doppler = fc*(max_speed*1000/3600)/c;
path_Doppler = max_Doppler*cos(2*pi*rand(1,taps)-pi); % Jakes' formula

rayleighchan = comm.RayleighChannel(...    
    'SampleRate', fs,...
    'PathDelays', delays,...
    'AveragePathGains', pdp,...
    'MaximumDopplerShift', 10E-20,...
    'PathGainsOutputPort', true...
    );
[~,pathgains] = rayleighchan(Tx_waveform);
for ii = 1:taps   %手動為各條path加入Doppler
    pathgains(:,ii) = pathgains(1,ii) * exp(1j*2*pi*((0:size(pathgains,1)-1).'/fs)*path_Doppler(ii));
end
chaninfo = info(rayleighchan);
coeff = chaninfo.ChannelFilterCoefficients;
delaydata = zeros(size(Tx_waveform,1),taps);
for ii = 1:taps
    delaydata(:,ii) = filter(coeff(ii,:),1,Tx_waveform);
end
Rx_waveform = sum(pathgains .* delaydata,2);   %手動計算訊號通過channel的結果
Rx_waveform = awgn(Rx_waveform, SNR, 'measured');
