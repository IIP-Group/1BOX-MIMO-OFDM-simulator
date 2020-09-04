% =========================================================================
% -- pilot sequence generator for 1-bit channel estimation
% -- APRIL 2019 (c) Seyed Hadi Mirfarshbafan (sm2675@cornell.edu)
% =========================================================================

function [XF, XT] = generate_pilot_sequence(par)


XT = zeros(par.U, par.Ntones + par.CPLength, par.T);
XF0 = zeros(par.U, par.Nu, par.T);
XF = zeros(par.U, par.Ntones, par.T);
% Frequency domain pilot matrix
switch (par.pilot)
    case 'Hadamard'
        if ~strcmp(par.channel_freq_type, 'frequency-flat')
            error('Hadamard pilot only defined for frequency-flat (at least for now)!')
        end
        X0 = hadamard(par.T);
        XF0 = sqrt(par.Es)*X0(1:par.U,:);      % this is the U*T training matrix (works only for frequency-flat)
    case 'Gaussian'
        XF0 = sqrt(par.Ntones/par.Nu)*sqrt(par.Es/2)*(randn(par.U, par.Nu, par.T)+1i*randn(par.U, par.Nu, par.T));
    case 'randQPSK'
        XF0 = sqrt(par.Ntones/par.Nu)*2*sqrt(par.Es/2)*((randi(2,par.U, par.Nu, par.T)-1.5) + 1i*(randi(2,par.U, par.Nu, par.T)-1.5));
    case 'diagonal-const'
        X0 = zeros(par.U, par.T);
        for t = 1:par.T
            X0(mod(t,par.U)+par.U*(mod(t,par.U)==0),t) = 2*sqrt(par.Es/2)*sqrt(par.U)*((randi(2)-1.5) + 1i*(randi(2)-1.5));
        end
        for w = 1:par.Nu
            XF0(:,w,:) = X0;
        end
    case 'diagonal'
        for t = 1:par.T
            XF0(mod(t,par.U)+par.U*(mod(t,par.U)==0),:,t) = 2*sqrt(par.Es/2)*sqrt(par.U)*sqrt(par.Ntones/par.Nu)*((randi(2,par.Nu,1)-1.5) + 1i*(randi(2,par.Nu,1)-1.5));
        end
    otherwise
        error('par.pilot not defiend!')
end

for t = 1:par.T
    
    XF(:,par.Ng/2+1:par.Nu+par.Ng/2,t) = XF0(:,:,t);
    % Convert to time domain
    XT(:,par.CPLength + 1:end,t) = sqrt(par.Ntones)*ifft(fftshift(squeeze(XF(:,:,t)).', 1)).';
    % Add cyclic prefix
    XT(:,1:par.CPLength,t) = XT(:,end - par.CPLength + 1:end, t);
end

if strcmp(par.channel_freq_type, 'FF')
    XT = XF0;
    XT = squeeze(XT);
end

end
