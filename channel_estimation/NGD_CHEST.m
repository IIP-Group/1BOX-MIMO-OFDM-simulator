% =========================================================================
% -- APRIL 2019 (c) Seyed Hadi Mirfarshbafan (sm2675@cornell.edu)
% -----------------------------------------------------
% Normalized Gradient Descent Channel Estimator for 1-bit OFDM systems
% =========================================================================
% XF must be the U*W*T frequency-domain training matrix
% Yq must be the B*W*T time-domain (after removing CP) training signals received at the B BS antennas
function HhatN = NGD_CHEST(par, XF, Yq, N0)

HhatN = zeros(par.B,par.U,par.Ntones);
t_neg = -5;
t_pos = 3;

gamma = sqrt(par.Es);
% Set the step size
kappa = par.kappa_CHEST;
epsilon = 0.001;    % Termination Threshold
% --- sigma
N0F = N0;
if par.UseFakeSNRCHEST
    if N0 < par.fakeN0_CHEST
        N0F = par.fakeN0_CHEST;
    end
end

sigma = sqrt(N0F);

% --------- TDMLE denoising variables -------------
% according to the following paper:
% "OFDM Channel Estimation Algorithm and ASIC Implementation:
% ************ ACHTUNG!:
% The DFT matrix F here is NOT normalized, to comply with the TDMLE paper.
P = par.Ng/2+1:par.Ng/2+par.Nu;
F = fft(eye(par.Ntones));
FL = F(:,1:par.ChannelTaps);
FP = FL(P,:);

Gradient = zeros(par.U, par.Ntones);
Z = zeros(par.T, par.Ntones);
%Hhat = 1/(gamma*sqrt(2)/sigma)*ones(par.B,par.U,par.Ntones);
%Hhat = zeros(par.B,par.U,par.Ntones);
Hhat = ZF_CHEST(par, XF, Yq);


    
for b = 1:par.B
    iteration = 0;
    signReal = real(squeeze(Yq(b,:,:)).');
    signImag = imag(squeeze(Yq(b,:,:)).');
    H_k = squeeze(Hhat(b,:,:));
    H_km1 = H_k;
    while ( norm(H_k(:) - H_km1(:)) >= epsilon*norm(H_km1(:)) || iteration == 0 ) && iteration < par.max_iterations_CHEST

        H_km1 = H_k;
        
        for tone = 1:par.Ntones
            XFs = squeeze(XF(:,tone,:)).';
            Z(:,tone) = XFs*H_k(:,tone);
        end
        alpha = gamma*sqrt(2*par.Ntones)/sigma*ifft(fftshift(Z.',1)).';
        
        alpha = signReal.*real(alpha) + 1i*signImag.*imag(alpha);
        gradient = omega_complex(alpha, t_neg, t_pos);
        gradient = gamma*(signReal.*real(gradient) + 1i*signImag.*imag(gradient));
        Gf = 1/sqrt(par.Ntones)*fftshift(fft(gradient.'),1).';
        for tone = 1:par.Ntones
            XFs = squeeze(XF(:,tone,:)).';
            Gradient(:,tone) = XFs'*Gf(:,tone);
        end
        H_k = H_k + kappa*Gradient;

        iteration = iteration + 1;
    end % end iteration loop
    
    % Apply TDMLE denoising
    if par.denoise_channel
        for u = 1:par.U
            H_k(u,par.Ng/2+1:par.Ng/2+par.Nu) = (FP*((FP'*FP)\(FP'*H_k(u,par.Ng/2+1:par.Ng/2+par.Nu).'))).';
        end
    end
    Hhat(b,:,:) = H_k;
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


end % end function

%% complex omega function
% NOTE: x should be a complex vector
function g = omega_complex(x, t_neg, t_pos)

    g = omega_real(real(x), t_neg, t_pos) + 1i*omega_real(imag(x), t_neg, t_pos);
end


%% real omega function
% NOTE: x must be a real vector
function g = omega_real(x, t_neg, t_pos)
    g = 1/sqrt(2*pi)*(exp(-x.^2/2)./normcdf(x));
    g(x > t_pos) = 0;
    g(x < t_neg) = -x(x < t_neg);
end