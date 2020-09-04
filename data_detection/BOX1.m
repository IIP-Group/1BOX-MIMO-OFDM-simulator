% =========================================================================
% Title       : 1BOX detector for massive-MU-MIMO-OFDM
% File        : BOX1.m
% -------------------------------------------------------------------------
% Description :
% Argument description:
% par: The structure holding the simulation parameters.
% Hf: The frequency-domain channel matrices.
% yq: 1-bit quantized received data. This shall be a matrix containing the
% data of all BS antennas over W time samples.
% N0: noise variance
% -------------------------------------------------------------------------
% Revision: 1
% Date: 2/9/2018
% -------------------------------------------------------------------------
% Author: Seyed Hadi Mirfarshbafan
% =========================================================================

function [iteration,shat,idxnML,bitnML] = BOX1(par,Hf,yq, N0)
   
    t_neg = -4;
    t_pos = 4;
    % Set the step size
    kappa = par.kappa;
    epsilon = 0.001;    % Termination Threshold
    % --- sigma
    N0F = N0;
    if par.UseFakeSNRDET
        if N0 < par.fakeN0_DET
            N0F = par.fakeN0_DET;
        end
    end
    sigma = sqrt(N0F);
    
    
    % Normalize channel estimates!
    for w = (1+par.Ng/2):(par.Nu + par.Ng/2)
        for u = 1:par.U
            hhat = Hf(:,u,w);
            Hf(:,u,w) = sqrt(par.Eh)*sqrt(par.B)*hhat/norm(hhat);
        end
    end
    
    Hf = fftshift(Hf, 3);   % this is to avoid fftshift operations inside the algorithm loop.
                            % this is more efficient and closer to hardware
                            % implmentation in which there is no fftshift
                            % operation in the algorithm. See
                            % MPGD_OFDM_fftshift.m for a version with
                            % fftshift operations inside the algorithm if
                            % that's more intuitive.
    
    shat = zeros(par.U, par.Nu);
    idxnML = zeros(par.U, par.Nu);
    bitnML = zeros(par.U, par.Q, par.Nu);
    Z = zeros(par.B, par.Ntones);
    Gradient = zeros(par.U, par.Ntones);
    
    % initial point for the vectors of users' data
    S0 = zeros(par.U, par.Ntones);
    S_k = S0;
    S_km1 = S0;
    iteration = 0;

    while ( norm(S_k(:) - S_km1(:)) >= norm(S_km1(:))*epsilon || iteration == 0 ) && iteration < par.max_iterations
        S_km1 = S_k;
        for tone = 1:par.Ntones
            Z(:,tone) = Hf(:,:,tone)*S_k(:,tone);
        end
        
        argphi = sqrt(2*par.Ntones)/sigma*ifft(Z.').';
        argphi = real(argphi).*real(yq) + 1i*imag(argphi).*imag(yq);
        Gt = omega_complex(argphi, t_neg, t_pos);
        Gt = real(yq).*real(Gt) + 1i*imag(yq).*imag(Gt);
        Gf = 1/sqrt(par.Ntones)*fft(Gt.').';
        for tone = 1:par.Ntones
            Gradient(:,tone) = Hf(:,:,tone)'*Gf(:,tone);
        end
        S_k = S_k + kappa*Gradient;
        
        % project the result of k-th iteration
        S_k = Proj(S_k, par);
        % Nullify the unused subcarriers
        S_k(:,par.Nu/2 + 1:par.Nu/2 + par.Ng) = 0;
        iteration = iteration + 1;
    end % end iteration loop
    
    
    S_k = fftshift(S_k, 2);  % this is to avoid fftshift operations inside the algorithm loop.
                            % this is more efficient and closer to hardware
                            % implmentation in which there is no fftshift
                            % operation in the algorithm. See
                            % MPGD_OFDM_fftshift.m for a version with
                            % fftshift operations inside the algorithm if
                            % that's more intuitive.
                            
    % -- normalize results
    SN = sqrt(par.U)*sqrt(par.Nu)*S_k/norm(S_k, 'fro');
    % -- compute outputs - Hard decision
    for tone = 1:par.Nu
        shat(:,tone) = SN(:,tone + par.Ng/2);
        [~,idxnML(:,tone)] = min(abs(shat(:,tone)*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
        bitnML(:,:,tone) = par.bits(idxnML(:,tone),:);
    end
        
end

%% Orthogonal projection
function xpr = Proj(x, par)
    Tm = max(abs(real(par.symbols)));
    xpr = sign(real(x)).*min(abs(real(x)),Tm) + 1i*sign(imag(x)).*min(abs(imag(x)),Tm);
    
end

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