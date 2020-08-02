function [BER,ffPower,fbPower] = ModuloSKenv(nbits,N,R,snrdB,DsnrdB,Petarget)
tic
PAMSize = R*N;
inbits = 1*(rand(PAMSize*ceil(nbits/PAMSize),1)>0.5);
Theta = PAMmodulate(inbits,PAMSize);
SNR  = 10^(snrdB/10);
Dsnr = 10^(DsnrdB/10);
if (N==1) % just one round - simple PAM without feeback
    X = Theta;
    Y = AWGNchan(X,snrdB);
    outbits = PAMdemodulate(Y,PAMSize);
else % more than one round, use feedback
    if isinf(DsnrdB) % noiseleff feedback. use classic SK
        [outbits,ffPower,fbPower] = SKscheme(Theta,PAMSize,N,SNR);
    else
        [outbits,ffPower,fbPower] = MSKscheme(Theta,PAMSize,N,SNR,Dsnr,Petarget);
    end
end

BER = mean(1*(outbits~=inbits));
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outbits,ffPower,fbPower] = SKscheme(Theta,PAMSize,N,SNR)
dispLog = 1;
X = Theta;
Y = AWGNchan(X,SNR);
ThetaHat = Y;
sigma2chan = 1/SNR;
sigma2n = 1/SNR;
alphan = 1/sqrt(sigma2n);
ffPower = zeros(1,N);
fbPower = zeros(1,N);
%printLog(dispLog,['Running SK iteration #',num2str(1)]);
ffPower(1) = mean(X.^2);
for ii = 2:N
%    printLog(dispLog,['Running SK iteration #',num2str(ii)]);
    ffPower(ii) = mean(X.^2);
    errn = ThetaHat-Theta;
    X = alphan*errn;
    Y = AWGNchan(X,SNR);
    % advance
    betan = sqrt(sigma2n/sigma2chan)*sqrt(SNR)/(1+SNR);
    ThetaHat = ThetaHat - Y*betan;
    sigma2n = sigma2n/(1+SNR);
    alphan = 1/sqrt(sigma2n);
end
outbits = PAMdemodulate(ThetaHat,PAMSize);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outbits,ffPower,fbPower] = MSKscheme(Theta,PAMSize,N,SNR,Dsnr,Petarget)

% the constants
dispLog = 0;
d = sqrt(12);
pm = Petarget/2/(N-1);
L = 1/3*(qfuncinv(pm/2))^2;
Psi3 = 1+3/2*L*pm*(N-1-2/N);
snrtag = SNR/Psi3;
Psi1 = 1+Psi3*L/Dsnr;
Fsnr = snrtag*Dsnr; % feedback SNR
Psi2 = 1/(1-L/Fsnr);
snrn = snrtag*(1+snrtag/Psi1/Psi2).^((1:N)-1);
sigma2n = 1./snrn;
alpha0 = sqrt(L/Psi3);
sigma2f = 1/snrtag; % feedforward
sigma2b = 1/Fsnr; % feedback
betan = sqrt(sigma2n/sigma2f)*sqrt(snrtag*(1-L/Fsnr))/(1+snrtag);
gamman = sqrt(1/L-sigma2b).*sqrt(snrn);
% first round is like SK
X = Theta;
Y = AWGNchan(X,SNR);
ThetaHat = Y;
totffPower = mean(X.^2);
totfbPower = 0;
ffPower = zeros(1,N);
fbPower = zeros(1,N);

printLog(dispLog,['Running Modulo-SK. Feedforward Iteration #',num2str(1)]);
ffPower(1) = mean(X.^2);
for ii = 2:N
    % the feedback channel
    Xb = ModD(gamman(ii-1)*ThetaHat,d);
    totfbPower = totfbPower+mean(Xb.^2);
    printLog(dispLog,['Running Modulo-SK. Feedback Iteration #',num2str(ii-1)]);
    fbPower(ii-1) = mean(Xb.^2);
    Yb = AWGNchan(Xb,Fsnr);    
    % process Terminal A
    errnA = 1/gamman(ii-1)*ModD(Yb-gamman(ii-1)*Theta,d);
    X = alpha0*gamman(ii-1)*errnA;
    totffPower = totffPower+mean(X.^2);
    printLog(dispLog,['Running Modulo-SK. Feedforward Iteration #',num2str(ii)]);
    ffPower(ii) = mean(X.^2);
    Y = AWGNchan(X,SNR);
    % process Terminal B
    errnB = betan(ii-1)*Y;
    ThetaHat = ThetaHat-errnB;
end
outbits = PAMdemodulate(ThetaHat,PAMSize);
%disp(['Average FF power =',num2str(totffPower/N),', Average FB power =',num2str(totfbPower/N)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = AWGNchan(X,snr)
Y = X + randn(size(X))/sqrt(snr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = PAMmodulate(inbits,PAMSize)
eta = sqrt(3/(2^(2*PAMSize)-1));
bitsyms = reshape(inbits,PAMSize,[]);
X = eta*(2.^(PAMSize:-1:1)*bitsyms - (2^PAMSize-1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outbits = PAMdemodulate(Y,PAMSize)
eta = sqrt(3/(2^(2*PAMSize)-1));
Xh = round((Y/eta + (2^PAMSize-1))/2);
Xh = min(2^PAMSize-1,max(0,Xh));
bitsmat = (dec2bin(Xh,PAMSize) - '0')';
outbits = bitsmat(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function md = ModD(x,d)
md = x-d*round(x/d);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printLog(dispLog,str)
if dispLog
    disp(str);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%