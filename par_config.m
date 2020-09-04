% -----------------------------------------------------
% -- Simple MIMO simulator with estimated LoS channel
% -- 2018 (c) studer@cornell.edu, Seyed Hadi Mirfarshbafan (sm2675@cornell.edu)
% -----------------------------------------------------
function par = par_config(setup)

fullpath = mfilename('fullpath');
simulator_path = strrep(fullpath, mfilename(), '');
addpath(genpath(simulator_path));
par.simulator_path = simulator_path;

%% Parameters

% simulation parameters according to the four plots in figure 3 of the
% paper.

switch (setup)
    case 'a'
        par.runId = 1; % simulation ID (used to reproduce results)
        par.B = 128; % number of BS antennas
        par.U = 8; % number of UEs
        par.mod = '8PSK'; % modulation scheme: 'BPSK','QPSK', '8PSK', '16QAM'
        par.SNRdB_list_L = {[-5 0 5 6], [-5 0 5 10 15], [-5 0 5 10 15 ], [-5 0 5 10 15], [-5 0 5 10 15]}; % list of SNR [dB] values to be simulated
    case 'b'
        par.runId = 2; 
        par.B = 128; 
        par.U = 8; 
        par.mod = '16QAM'; 
        par.SNRdB_list_L = {[-5 0 5 8], [-5 0 5 10 15], [-5 0 5 10 15 ], [-5 0 5 10 15], [-5 0 5 10 15]}; 
    case 'c'
        par.runId = 10; 
        par.B = 64; 
        par.U = 4; 
        par.mod = '8PSK';
        par.SNRdB_list_L = {[-5 0 5 6], [-5 0 5 10 15], [-5 0 5 10 15 ], [-5 0 5 10 15], [-5 0 5 10 15]}; 
    case 'd'
        par.runId = 11; 
        par.B = 64; 
        par.U = 4; 
        par.mod = '16QAM'; 
        par.SNRdB_list_L = {[-5 0 5 8], [-5 0 5 10 15], [-5 0 5 10 15 ], [-5 0 5 10 15], [-5 0 5 10 15]}; 
    otherwise
        error('Input to the par_config function is not defined. It must be ''a'', or ''b'' or ''c'' or ''d''.');
end

par.quant_list = {inf, 1,1,1,1};
par.chest_list = {'ZF-CHEST', 'ZF-CHEST', 'PerfectCSI', 'ZF-CHEST', 'NGD-CHEST'}; 
par.det_list = {'ZF-DET', 'ZF-DET', '1BOX', '1BOX', '1BOX'};


if length(par.quant_list) ~= length(par.chest_list) || length(par.quant_list) ~= length(par.det_list)
    error('number of elements of the cell arrays par.quant_list and par.chest_list and par.det_list, must be the same.')
end


par.channel_type = 'Rayleigh';
par.chest_normalization_method = 'per_BS_antenna'; 
par.denoise_channel = 1;
par.T = 2*par.U;    % channel training length, NOTE: we must have: par.T >= par.U
%WARNING: par.T must be a multiple of par.U
par.pilot = 'randQPSK';    

% -------------------- Detector parameters -----------------------------
par.kappa = 1/64*sqrt(2);               % For ease of implementation (1/128 + 1/256 is also good)
par.kappaFP = 1/32;
par.kappa_CHEST = 1/16;
par.fakeSNR = 10;   % dB
par.fakeSNRFP = 10;   % dB

par.UseFakeSNRDET = 1;          % Use fake SNR in 1BOX detector
par.UseFakeSNRDETFP = 1;        % Use fake SNR in fixed point detector
par.UseFakeSNRCHEST = 0;        % Use fake SNR in NGD channel estimator

par.max_iterations = 3;
par.max_iterations_CHEST = 5;

% -------------------- Simulation parameters -----------------------------


par.minVER = 100; % Minimum number of vector errors to see before ending trials for each SNR

% ######## not needed to modify:

par.channel_freq_type = 'FS';   % 'FS' or 'FF'
par.Eh = 1; % Variance of each entry of the channel matrix
par.fakeN0_CHEST = par.Eh*par.U*10^(-par.fakeSNR/10);
par.fakeN0_DET = par.Eh*par.U*10^(-par.fakeSNR/10);
par.fakeN0_DETFP = par.Eh*par.U*10^(-par.fakeSNRFP/10);
par.parallel_trials = 1e2; % number of Monte-Carlo trials within each while loop
par.max_trials = 5000;
% ----------------------- Parameters of OFDM --------------------------

if strcmp(par.channel_freq_type, 'FF')
    par.ChannelTaps = 1;
    par.Nu = 1;
    par.Ng = 0;
    par.Ntones = 1;
    par.CPLength = 0;
else
    par.ChannelTaps = 15;
    par.Nu = 100;           % Number of OFDM tones used for data
    par.Ng = 28;            % Number of OFDM guard tones
    par.Ntones = par.Nu + par.Ng;
    par.CPLength = par.ChannelTaps - 1;
end


%% Initializations : DO NOT MODIFY
% -------------------------------------------------------------------
% set up Gray-mapped constellation alphabet (according to IEEE 802.11)
switch (par.mod)
    case 'BPSK'
        par.symbols = [ -1 1 ];
    case 'QPSK'
        par.symbols = [ -1-1i,-1+1i, ...
            +1-1i,+1+1i ];
    case '8PSK'
        par.symbols = [sqrt(2), 1+1i, -1 + 1i, 1i*sqrt(2), 1 - 1i,...
            -1i*sqrt(2), -sqrt(2), -1-1i ];
    case '16QAM'
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    otherwise
        error('modulation not defined!');
        
end

% extract average symbol energy
par.Es = mean(abs(par.symbols).^2);
par.symbols = par.symbols/sqrt(par.Es);
par.Es = mean(abs(par.symbols).^2);

% precompute bit labels
par.Q = log2(length(par.symbols)); % number of bits per symbol
par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');

end