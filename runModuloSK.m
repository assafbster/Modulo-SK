% examples for using the feedback schemes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% General setting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
% the prefered SK setting
% DsnrdB = inf;    % the ratio between feedback SNR and feedforward SNR
%                  % inf means clean feedback and classical SK scheme
% N = 15;
% R = 1/3; 

% the prefered ModuloSK setting
DsnrdB = 20;      % the ratio between feedback SNR and feedforward SNR
N = 24;           % number of SK iterations   
R = 1/3;          % rate

Petarget = 1e-6; % the target BER
                 % Note that systemwise it is used as the target error in
                 % the decoding a PAM symbol. So BER is smaller for PAM
                 % modulations with more than one bit

nPAMsyms = 1e6;      % number of PAM symbols
nbits = nPAMsyms*R*N % number of simulated bits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find an SNR which provides the target BER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[snrShannondB,CapGapdB] = calcSNRworkPoint(N,R,DsnrdB,Petarget);
snrdB  = snrShannondB + CapGapdB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Typing announcements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['R = ',num2str(R),', Feedforward SNR =',num2str(snrdB),'dB, which is ',...
    num2str(CapGapdB),'dB more than Shannon''s limit']);
if isinf(DsnrdB)
    disp(['Noiseless Feedback. Using classical SK for ',num2str(N),...
        ' iterations']);
else
    disp(['Feedforward SNR =',num2str(snrdB),'. Using Modulo-SK for ',num2str(N),...
        ' iterations']);
end
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Running the coding scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BER = ModuloSKenv(nbits,N,R,snrdB,DsnrdB,Petarget); % PE is BER
disp(['Simulation ended. Total measured BER = ',num2str(BER)]);

