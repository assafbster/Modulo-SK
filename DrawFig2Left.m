starttime = tic;
nPAMsyms = 100e6;      % number of PAM symbols
nPAMbatch = 1e6;

R = 1/3;          % rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find an SNR which provides the target BER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nvec = [150,39,15];
bsnrvec = [inf,inf,inf];
snrdBVecCell =  cell(length(Nvec),1);
snrdBVecCell{1} = -2:0.1:-1.8;
snrdBVecCell{2} = -1.7:0.1:-1.2;
snrdBVecCell{3} = -0.8:0.1:0;
colorcell = {'k','b','r','c'};
legendcell = cell(2*length(Nvec),1);
for kk = 1:length(Nvec)
    bsnr = bsnrvec(kk);
    N = Nvec(kk); 
    legendcell{2*kk-1} = ['K = ',num2str(N/3),' predicted'];
    legendcell{2*kk} = ['K = ',num2str(N/3),' simulated'];
    snrdBVec = snrdBVecCell{kk};
    DsnrdBVec = bsnr-snrdBVec;
    PeTargetPerDSNR = zeros(size(DsnrdBVec));
    maxerr = qfunc(sqrt(3));
    PeTargetVec = logspace(-3,-20,1000);
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
        [~,maxind] = max(find(chanSNRvec<snrdBVec(ii)));
        if isempty(maxind)
            PeTargetPerDSNR(ii) = maxerr;
            disp('no success');
        else
            PeTargetPerDSNR(ii) = PeTargetVec(maxind);
        end
    end
    
    semilogy(snrdBVec,PeTargetPerDSNR,[colorcell{kk},'--']);
    hold on;
    grid on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BERvec = zeros(size(DsnrdBVec));
    nbatches = ceil(nPAMsyms/nPAMbatch);
    curBERvec = zeros(nbatches,1);
    for ii = 1:length(DsnrdBVec)
        disp('####################################################################');
        disp(['snrdBVec = ',num2str(snrdBVec(ii))]);
        disp('####################################################################');
        curBERvec = zeros(1,nbatches);
        ffPowerMat = zeros(nbatches,N);
        fbPowerMat = zeros(nbatches,N);
        for jj = 1:nbatches
            disp(['batch # = ',num2str(jj),' out of ',num2str(nbatches)]);
            nbits = nPAMbatch*R*N; % number of simulated bits
            [curBERvec(jj),ffPowerMat(jj,:),fbPowerMat(jj,:)] = ModuloSKenv(nbits,N,R,snrdBVec(ii),DsnrdBVec(ii),PeTargetPerDSNR(ii));
        end
        disp(['curBERvec = ',num2str(curBERvec)]);
        disp(['ff power per iter = ',num2str(mean(ffPowerMat,1))]);
%         disp(['fb power per iter = ',num2str(mean(fbPowerMat,1))]);
        disp(['average ff power = ',num2str(mean(mean(ffPowerMat)))]);
%         disp(['average fb power = ',num2str(mean(mean(fbPowerMat)))]);
        BERvec(ii) = mean(curBERvec);
    end
    semilogy(snrdBVec,BERvec,[colorcell{kk},'-']);

end
axis([-2,1,1e-20,1e-3])
legend(legendcell,'Location','northeast');
xlabel('SNR on feedforward channel','FontSize',24);
ylabel('BER','FontSize',24);

toc(starttime)
% save('DrawFig2Left.mat');