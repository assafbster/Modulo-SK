function [snrShannondB,CapGapdB,pointFound] = calcSNRworkPoint(N,R,DsnrdB,Petarget)

Dsnr = 10^(DsnrdB/10);

PAMcapGap = 1/3*(qfuncinv(Petarget/2))^2;
PAMcapGapdB = 10*log10(PAMcapGap);
snrShannondB = 10*log10(2^(2*R)-1);
snrdBmin = snrShannondB+min(0,PAMcapGapdB);
snrdBmax = snrShannondB+max(9,PAMcapGapdB);
if (N==1)
    CapGapdB = PAMcapGapdB;
else
    if  isinf(DsnrdB) % noiseleff feedback. use classic SK
        Gamma0SK = 1/3*(qfuncinv(Petarget/2))^2;
        snrVecdB = snrdBmin:0.01:snrdBmax;
        for ii = 1:length(snrVecdB)
            SNR = 10^(snrVecdB(ii)/10);
            snrn = SNR*(1+SNR)^(N-1);
            achRate = 0.5/N*log2(1+snrn/Gamma0SK);
            if achRate>R
                pointFound = 1;
                break
            end
        end
        CapGapdB = 10*log10(SNR)-snrShannondB;        
    else % Modulo SK
        Gamma0 = 1/3*(qfuncinv(Petarget/4))^2;
        pm = Petarget/2/(N-1);
        L = 1/3*(qfuncinv(pm/2))^2;
        Psi3 = 1+3/2*L*pm*(N-1-2/N);
        Psi1 = 1+Psi3*L/Dsnr;
        snrVecdB = snrdBmin:0.01:snrdBmax;
        pointFound = 0;
        for ii = 1:length(snrVecdB)
            SNR = 10^(snrVecdB(ii)/10);
            snrtag = SNR/Psi3;
            Fsnr = snrtag*Dsnr; % feedback SNR
            Psi2 = 1/(1-L/Fsnr);
            if Psi2>0
                snrn = snrtag*(1+snrtag/Psi1/Psi2)^(N-1);
                achRate = 0.5/N*log2(1+snrn/Gamma0);
                if achRate>R
                    pointFound = 1;
                    break
                end
            end
        end
        if pointFound
            CapGapdB = 10*log10(SNR)-snrShannondB;
        else
            CapGapdB = PAMcapGapdB;
        end
    end    
end