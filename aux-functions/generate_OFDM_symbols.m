% =========================================================================
% Title       : Simulator for Quanized Massive MU-MIMO-OFDM Uplink
% File        : generate_data
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

function [bits, idx, S, xt] = generate_OFDM_symbols(par)
  
  % Generate bits for all users and all subcarriers
  bits = randi([0 1],par.U,par.Q, par.Nu);
  idx = ones(par.U, par.Nu);
  % generate transmit symbol
  S = zeros(par.Ntones, par.U);
  x = zeros(par.Ntones + par.CPLength, par.U);
  for ii = 1:par.Nu
    idx(:,ii) = bi2de(bits(:,:,ii),'left-msb')+1;
  end
  S(par.Ng/2+1:par.Nu+par.Ng/2, :) = sqrt(par.Ntones/par.Nu)*par.symbols(idx).';   % This is the frequency domain symbol set of all users
                            % Each column of S, contains the frequency domain
                            % symbols of the corresponding user.
  % Convert to time domain
  x(par.CPLength + 1:end,:) = sqrt(par.Ntones)*ifft(fftshift(S, 1));   
  % Add cyclic prefix
  x(1:par.CPLength,:) = x(end - par.CPLength + 1:end, :);
  
  if par.ChannelTaps == 1
      x = S;
  end
  
  xt = x.';
  
end
