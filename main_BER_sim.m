% -----------------------------------------------------
% -- main simulator for 1-bit detection in massive MU-MIMO-OFDM systems
% -- MARCH 2019 (c) Seyed Hadi Mirfarshbafan (sm2675@cornell.edu)
% -----------------------------------------------------
clearvars;

simulation_setup = 'd'; % To regenerate the plots in figure 3 of the paper use one of the following: 'a', 'b', 'c', 'd'.
% To simulate with other parameters, define set of parameters in
% par_config.m and then initialize simulation_setup accordingly.

par = par_config(simulation_setup);

% use runId random seed (enables reproducibility)
rng(par.runId, 'twister');
tic
% -- start simulation



Hf_S = zeros(par.B,par.U,par.Ntones, par.parallel_trials);
H_S = zeros(par.B,par.U,par.Ntones, par.parallel_trials);
bits_S = zeros(par.U, par.Q, par.Nu, par.parallel_trials);
idx_S = zeros(par.U, par.Nu, par.parallel_trials);
S_S = zeros(par.parallel_trials,par.Ntones, par.U);
x_S = zeros(par.U, par.Ntones + par.CPLength, par.parallel_trials);
N_CHEST_S = zeros(par.B,par.Ntones, par.T,par.parallel_trials);
N_S = zeros(par.B,par.Ntones,par.parallel_trials);

[XF, XT] = generate_pilot_sequence(par);

%% SNR loop

for sim_idx = 1:length(par.det_list)
    
    par.detector = par.det_list{sim_idx};
    par.ch_estimator = par.chest_list{sim_idx};
    par.b = par.quant_list{sim_idx};
    par.SNRdB_list = par.SNRdB_list_L{sim_idx};
    
    % initialize result arrays (SNR x 1)
    res.VER = zeros(length(par.SNRdB_list),1); % vector error rate
    res.SER = zeros(length(par.SNRdB_list),1); % symbol error rate
    res.BER = zeros(length(par.SNRdB_list),1); % bit error rate
    res.min_trials = zeros(length(par.SNRdB_list),1);
    
    VERtemp = zeros(par.parallel_trials,1);
    SERtemp = zeros(par.parallel_trials,1);
    BERtemp = zeros(par.parallel_trials,1);
    
    fprintf(['Starting simulation for B = ' num2str(par.B) ', U = ' num2str(par.U) ', ' num2str(par.b) ' bit ADCs, channel estimator: ' par.ch_estimator ', detector: ' par.detector ' ...\n']);
    
    for snr_idx = 1:length(par.SNRdB_list)
        
        % compute noise variance (average SNR per receive antenna is: SNR=U*Es*Eh/N0)
        N0 = par.Eh*par.U*par.Es*10^(-par.SNRdB_list(snr_idx)/10);
        % trials loop: repeat until enough vector errors occur
        while (res.VER(snr_idx) < par.minVER) && (res.min_trials(snr_idx)*par.parallel_trials < par.max_trials)
            
            % Create the channel matrix, users' data to be transmitted and
            % noise matrices. This is done seperately from the parfor loop to
            % enable fixed random seed and reprodocibility.
            for t_idx=1:par.parallel_trials
                [Hf_S(:,:,:,t_idx), H_S(:,:,:,t_idx)] = generate_channels(par);     % Hf is frequency domain and H is time domain
                
                % generate one OFDM symbol per user. x is a matrix, each row
                % of which contains the OFDM time domain symbol for user
                % corresponding to that row.
                [bits_S(:,:,:,t_idx), idx_S(:,:,t_idx), S_S(t_idx,:,:), x_S(:,:,t_idx)] = generate_OFDM_symbols(par);
                
                % generate noise matrix at the receive antenna array, for the
                % channel estimation phase
                N_CHEST_S(:,:,:,t_idx) = sqrt(N0/2)*(randn(par.B, par.Ntones, par.T)+1i*randn(par.B, par.Ntones, par.T));
                
                % generate noise matrix at the receive antenna array
                N_S(:,:,t_idx) = sqrt(N0/2)*(randn(par.B,par.Ntones)+1i*randn(par.B,par.Ntones));
            end
            
            %  Parallel for loop
            
            
            
            parfor t_idx=1:par.parallel_trials
                % Create the channel matrix based on the chosen model
                Hf = Hf_S(:,:,:,t_idx);
                H = H_S(:,:,:,t_idx);
                
                
                bits = bits_S(:,:,:,t_idx);
                idx = idx_S(:,:,t_idx);
                S = S_S(t_idx,:,:);
                x = x_S(:,:,t_idx);
                
                % ---------------- channel training phase ----------------
                % generate noise matrix at the receive antenna array, for the
                % channel estimation phase
                N_CHEST = N_CHEST_S(:,:,:,t_idx);
                
                % transmit pilots over noisy channel, receive the unquantized
                % signal at the BS and remove CP.
                if ~strcmp(par.ch_estimator, 'PerfectCSI')
                    Y_train = zeros(par.B, par.Ntones, par.T);
                    for n_train = 1:par.T
                        
                        for nsym = 1:par.Ntones
                            for tap = 1:par.ChannelTaps
                                Y_train(:, nsym, n_train) = Y_train(:, nsym, n_train) + H(:,:,tap)*XT(:,nsym + par.ChannelTaps - tap, n_train);
                            end
                        end
                        
                    end
                    
                    % IMPORTANT: Here we assume that
                    % E[trace(XF*XF')] = par.Es*par.T*par.U
                    Y_train = 1/sqrt(par.Ntones)*Y_train;
                    Y_train = Y_train + N_CHEST;
                    
                    if par.b == 1
                        Yq_train = quantize(Y_train);
                    else
                        Yq_train = Y_train;
                    end
                end
                switch (par.ch_estimator)
                    case 'PerfectCSI'
                        Hest = Hf;
                    case 'Vanilla'
                        Hest = quantize(sqrt(par.Es*par.U)*Hf + sqrt(N0/2)*(randn(size(Hf))+1i*randn(size(Hf))));
                    case 'NGD-CHEST'
                        Hest = NGD_CHEST(par, XF, Yq_train, N0);
                    case 'ZF-CHEST'
                        Hest = ZF_CHEST(par, XF, Yq_train);
                    otherwise
                        error('par.denoiser not defiend')
                end
                
                
                % ---------------- data transmission phase ----------------
                
                N = N_S(:,:,t_idx);
                % transmit data over noisy channel, receive the unquantized
                % signal at the BS and remove CP.
                y = zeros(par.B, par.Ntones);
                for nsym = 1:par.Ntones
                    for tap = 1:par.ChannelTaps
                        y(:, nsym) = y(:, nsym) + H(:,:,tap)*x(:,nsym + par.ChannelTaps - tap);
                    end
                end
                y = 1/sqrt(par.Ntones)*y;
                y = y + N;
                % quantization at the receiver
                if par.b == 1
                    yq = quantize(y);
                else
                    yq = y;
                end
                
                switch (par.detector)     % select algorithms
                    
                    case 'ZF-DET'
                        [xhatN,idxhat,bithat] = ZF_DET(par,Hest,yq);
                    case '1BOX'
                        [iterations,xhat,idxhat,bithat] = BOX1(par,Hest,yq, N0);
                        
                    otherwise
                        error('par.detector type not defined.')
                end
                
                % -- compute error metrics
                err = (idx~=idxhat);
                VERtemp(t_idx) = sum(any(err)); % it should be scaled by 1/par.Nu for correct results
                SERtemp(t_idx) = sum(sum(err))/(par.U*par.Nu);
                BERtemp(t_idx) = sum(sum(sum(bits ~= bithat)))/(par.U*par.Nu*par.Q);
                
                
            end % end parfor
            res.VER(snr_idx) = res.VER(snr_idx) + sum(VERtemp);
            res.SER(snr_idx) = res.SER(snr_idx) + sum(SERtemp);
            res.BER(snr_idx) = res.BER(snr_idx) + sum(BERtemp);
            res.min_trials(snr_idx) = res.min_trials(snr_idx) + 1;
            fprintf('SNR: %d dB\tnumber of trials performed: %d\n', par.SNRdB_list(snr_idx), res.min_trials(snr_idx)*par.parallel_trials);
        end % end while trials loop
        
    end % SNR loop
    sim_name{sim_idx} = ['BER-' num2str(par.runId) '_' num2str(par.B) 'x' num2str(par.U) '_' par.mod '_Q' num2str(par.b) par.ch_estimator '_' par.detector]; % simulation name (used for saving results)
    legend_name{sim_idx} = [par.ch_estimator '/' par.detector];
    if isinf(par.b)
        legend_name{sim_idx} = [legend_name{sim_idx} ' (inf. res.)'];
    end
    save([par.simulator_path 'results/sim_res/' sim_name{sim_idx}],'par','res');
end % end simulation loop

% -- save final results (par and res structure)

plot_BER(1,par.runId,sim_name, legend_name)
toc