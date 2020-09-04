% =========================================================================
% Title       : Simulator for Quanized Massive MU-MIMO-OFDM Uplink
% File        : channel.m
% -------------------------------------------------------------------------
% Description :
% Generates the time domain and frequency domain channels using parameters
% given in par.
% -------------------------------------------------------------------------
% Revision: 0
% Date: 1/21/2018
% -------------------------------------------------------------------------
% Author: Seyed Hadi Mirfarshbafan
% =========================================================================

function [Hf, Ht] = generate_channels(par)

Hf = zeros(par.B,par.U,par.Ntones);
Ht = zeros(par.B,par.U,par.Ntones);

% -- generate channels in time domain
switch (par.channel_type)
    case 'Rayleigh'
        Ht(:,:,1:par.ChannelTaps) = sqrt(par.Eh)*sqrt(par.Ntones)*sqrt(0.5/par.ChannelTaps)*(randn(par.B,par.U,par.ChannelTaps) ...
            + 1i*randn(par.B,par.U,par.ChannelTaps));
     
    otherwise
        error('par.channel not defined YET!')
end




% -- convert to frequency domain
if strcmp(par.channel_freq_type, 'FF')
    Hf = Ht;
else
    for kk = 1:par.B
        for ll = 1:par.U
            Hf(kk,ll,:) = fftshift(fft(squeeze(Ht(kk,ll,:))))/sqrt(par.Ntones);
        end
    end
end
  
end
