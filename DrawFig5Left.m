starttime = tic;
nPAMsyms = 100e6;      % number of PAM symbols
nPAMbatch = 1e6;

N = 39;           % number of SK iterations
R = 1/3;          % rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find an SNR which provides the target BER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DsnrdBVec = [14.2,14.5:0.5:20];      % the ratio between feedback SNR and feedforward SNR

PeTargetPerDSNR = zeros(size(DsnrdBVec));
maxerr = qfunc(sqrt(3));
PeTargetVec = logspace(log10(maxerr),-9,1000);
% for Petarget
for ii = 1:length(DsnrdBVec)
    chanSNRvec = zeros(size(PeTargetVec));
    for jj = 1:length(PeTargetVec)
        [snrShannondB,CapGapdB,success] = calcSNRworkPoint(N,R,DsnrdBVec(ii),PeTargetVec(jj));
        if success
            chanSNRvec(jj) = snrShannondB + CapGapdB;
        else
            chanSNRvec(jj) = inf;
        end
    end
    [~,maxind] = max(find(chanSNRvec<0));
    if isempty(maxind)
        PeTargetPerDSNR(ii) = maxerr;
    else
        PeTargetPerDSNR(ii) = PeTargetVec(maxind);
    end
end

semilogy(DsnrdBVec,PeTargetPerDSNR);
hold on;
grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BERvec = zeros(size(DsnrdBVec));
snrdB = 0;
nbatches = ceil(nPAMsyms/nPAMbatch);
for ii = 1:length(DsnrdBVec)
    disp('##################################');
    disp(['DsnrdBVec = ',num2str(DsnrdBVec(ii))]);
    curBERvec = zeros(nbatches,1);
    for jj = 1:nbatches
    disp(['batch # = ',num2str(jj),' out of ',num2str(nbatches)]);
     nbits = nPAMbatch*R*N; % number of simulated bits
     curBERvec(jj) = ModuloSKenv(nbits,N,R,snrdB,DsnrdBVec(ii),PeTargetPerDSNR(ii)); 
    end
    BERvec(ii) = mean(curBERvec);
end
hold on;
semilogy(DsnrdBVec,BERvec,'r');
semilogy([0,20],[1,1]*1e-3/50,'k-.');
semilogy([0,20],[1,1]*0.0012,'k--');
axis([0,20,1e-7,1])
legend('Predicted BER','Simulated BER','metaconverse for K = 50','metaconverse for K = 13');
xlabel('SNR of feedback channel','FontSize',24);
ylabel('BER','FontSize',24);
toc(starttime)
save('DrawFig5Left.mat');
% nPAMsyms = 1e7 runtime is 346 seconds and BER of 1e-6 is smooth
% Take nPAMsyms = 1e8

% Getting the meta converse bound without feedback addpath spectre-master\awgn\
% converse(150,1e-3,1)/150 = 0.333 % for n = 150. the BER is lower bounded
% by 1e-3/50


% for n = 39
% ep = 0.0155;n = 39;converse(n,ep,1)/n
% ans = 0.3334
% ep/(n/3)=  0.0012



