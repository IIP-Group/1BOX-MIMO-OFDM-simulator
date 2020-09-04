%% ZF detector for OFDM
function [xhatN,idxhat,bithat] = ZF_DET(par,Hf,yq)
  yf = 1/sqrt(par.Ntones)*fftshift(fft(yq.'), 1).';
  if par.ChannelTaps == 1       % If single-tap then no need to convert the time domain signal back to frequency domain.
      yf = yq;                  % We'll work with the time domain signal instead, as was the case in the frequency-flat channel.
  end                           % Don't worry about the channel argument (Hf)! It is already kept in time domain.
  xhat = zeros(par.U, par.Nu);
  xhatN = zeros(par.U, par.Nu);
  idxhat = zeros(par.U, par.Nu);
  bithat = zeros(par.U, par.Q, par.Nu);
  for v = 1:par.Nu
      vv = v + par.Ng/2;
      Z = (Hf(:,:,vv)'*Hf(:,:,vv));
      W = Z\Hf(:,:,vv)';
      xhat(:,v) = W*yf(:,vv);
      
  end
  
  xhatN = xhat;
  if par.b
      xhatN = sqrt(par.U)*sqrt(par.Nu)*xhat/norm(xhat, 'fro');
  end
  
  for v = 1:par.Nu
      [~,idxhat(:,v)] = min(abs(xhatN(:,v)*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
      bithat(:,:,v) = par.bits(idxhat(:,v),:);
  end   
  
end