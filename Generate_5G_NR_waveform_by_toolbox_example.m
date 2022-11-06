clear; close all; clc;
% filter
numFFT = 1024;           % Number of FFT points
numRBs = 216;             % Number of resource blocks
rbSize = 12;             % Number of subcarriers per resource block
cpLen = 72;              % Cyclic prefix length in samples

bitsPerSubCarrier = 4;   % 2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM
snrdB = 20;              % SNR in dB

toneOffset = 2.5;        % Tone offset or excess bandwidth (in subcarriers)
L = 513;                 % Filter length (=filterOrder+1), odd

numDataCarriers = numRBs*rbSize;    % number of data subcarriers in sub-band
halfFilt = floor(L/2);
n = -halfFilt:halfFilt;

% Sinc function prototype filter
pb = sinc((numDataCarriers+2*toneOffset).*n./numFFT);

% Sinc truncation window
w = (0.5*(1+cos(2*pi.*n/(L-1)))).^0.6;

% Normalized lowpass filter coefficients
fnum = (pb.*w)/sum(pb.*w);

% Filter impulse response
% h = fvtool(fnum, 'Analysis', 'impulse', 'Fs', 15.36e6);

% Use dsp filter objects for filtering
filtTx = dsp.FIRFilter('Structure', 'Direct form symmetric', ...
    'Numerator', fnum);
filtRx = clone(filtTx); % Matched filter for the Rx

SNR = 20;
nrb = 216;
carrier = nrCarrierConfig('SubcarrierSpacing', 15, 'NSizeGrid', nrb);

% PDSCH
ncellid = 1007;
pdsch = nrPDSCHConfig;
pdsch.PRBSet = 0: carrier.NSizeGrid-1;      % Allocate the complete carrier
pdsch.SymbolAllocation = [0 14] ;           % Symbol allocation [Start Length]
pdsch.MappingType = 'A';                    % PDSCH mapping type ('A' or 'B')
pdsch.DMRS.DMRSTypeAPosition = 2;           % 2 or 3     
pdsch.DMRS.DMRSLength = 1;                  % 1 or 2
pdsch.DMRS.DMRSAdditionalPosition = 3;      % 0...3
pdsch.DMRS.DMRSConfigurationType = 1;       % 1 or 2
pdsch.DMRS.DMRSPortSet = 0;
pdsch.DMRS.NumCDMGroupsWithoutData = 1;     % 1 corresponds to CDM group number 0
pdsch.NumLayers = numel(pdsch.DMRS.DMRSPortSet);
pdsch.Modulation = '16QAM';
[ind,infor] = nrPDSCHIndices(carrier,pdsch,'IndexStyle','index','IndexOrientation','carrier');
numDataBits = infor.G;
cws = randi([0 1],numDataBits,1);
sym = nrSymbolModulate(cws,'16QAM');
sym_DMRS = nrPDSCHDMRS(carrier,pdsch, 'OutputDataType','double');
ind_DMRS = nrPDSCHDMRSIndices(carrier,pdsch,'IndexBase','1based','IndexOrientation','carrier');


% time-frequency grid
grid = nrResourceGrid(carrier);
grid(ind_DMRS) = sym_DMRS;
figure;
imagesc(abs(grid(1:12,:)));
xlabel('OFDM Symbols');
ylabel('Subcarriers');
title('DMRS allocation in one resource block');
set(gca,'FontSize',16)
grid(ind) = sym;

grid = fft(grid);

Tx_waveform = nrOFDMModulate(carrier, grid);

Tx_waveform = filtTx([Tx_waveform; zeros(L-1,1)]);


% channel
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


% receive
Rx_waveform = filtRx(Rx_waveform);
Rx_waveform = Rx_waveform(L:end);

grid_DM = nrOFDMDemodulate(carrier, Rx_waveform);
grid_DM = ifft(grid_DM);
grid_DM = fftshift(grid_DM);
[H, nVar, ~] = nrChannelEstimate(grid_DM,ind_DMRS,sym_DMRS);
[eqSym, ~] = nrEqualizeMMSE(reshape(grid_DM, [], 1), reshape(H, [], 1), nVar);
sym_Extract = nrExtractResources(ind, eqSym);

% figure;
% imagesc(abs(sym_Extract(1:12,:)));
% xlabel('OFDM Symbols');
% ylabel('Subcarriers');
% title('DMRS allocation in one resource block');
% set(gca,'FontSize',16)

demodbits = nrSymbolDemodulate(sym_Extract,'16QAM','DecisionType','Hard');
BER = sum(demodbits~=cws,'all') / numDataBits;

% PLOT
figure;                                                   
plot( reshape(sym_Extract, [], 1), '.' );
xlim([-2 2])
ylim([-2 2])
set(gca,'FontSize',16)
title(['signal constellation at Rx, BER = ', num2str(BER)]);

