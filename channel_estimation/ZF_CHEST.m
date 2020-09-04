%% ZF channel estimator for OFDM 
% XF must be the U*W*T frequency-domain training matrix
% Yq must be the B*W*T time-domain (after removing CP) training signals received at the B BS antennas
function HhatN = ZF_CHEST(par,XF,Yq)

HhatN = zeros(par.B,par.U,par.Ntones);
Hhat = zeros(par.B,par.U,par.Ntones);
Yqf = zeros(size(Yq));


% --------- TDMLE denoising variables -------------
% according to the following paper:
% "OFDM Channel Estimation Algorithm and ASIC Implementation:
% ************ ACHTUNG!:
% The DFT matrix F here is NOT normalized, to comply with the paper.
P = par.Ng/2+1:par.Ng/2+par.Nu;
F = fft(eye(par.Ntones));
FL = F(:,1:par.ChannelTaps);
FP = FL(P,:);

%% Take FFT
%######### IMPORTANT: ###########
% Note the 1/(par.Ntones) scale before fft! It is to make sure it is 
for b = 1:par.B
    Yqf(b,:,:) = 1/(par.Ntones)*fftshift(fft(squeeze(Yq(b,:,:))), 1);
end


%% Apply ZF filter per subcarrier
for w = (1+par.Ng/2):(par.Nu + par.Ng/2)
    XFs = squeeze(XF(:,w,:)).';
    W = (XFs'*XFs)\(XFs');
    for b = 1:par.B
        Yqfbw = squeeze(Yqf(b,w,:));
        Hhat(b,:,w) = W*Yqfbw;
    end
end


%% Apply TDMLE denoising
if par.denoise_channel  
    
    for b = 1:par.B
        for uu = 1:par.U
            hh = squeeze(Hhat(b, uu, par.Ng/2+1:par.Ng/2+par.Nu));
            Hhat(b, uu, par.Ng/2+1:par.Ng/2+par.Nu) = (FP*((FP'*FP)\(FP'*hh))).';
            
        end
    end
end
    

%% Normalize the results
if strcmp(par.chest_normalization_method, 'per_UE')
    for u = 1:par.U
        HhatN(:,u,:) = sqrt(par.Eh)*Hhat(:,u,:)/sqrt(mean(mean(abs(Hhat(:,u,:)).^2)));
    end
else
    for b = 1:par.B
        hhat = squeeze(Hhat(b,:,:));
        HhatN(b,:,:) = sqrt(par.Eh)*sqrt(par.U)*sqrt(par.Nu)*hhat/norm(hhat(:));
    end
end

  
end