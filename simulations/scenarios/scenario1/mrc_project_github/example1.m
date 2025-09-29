clear;clc;
%by simulation
rayChan = comm.RayleighChannel( ...
        'SampleRate',10000, ...
        'MaximumDopplerShift',100);
outage_values_SC=[];
outage_values_MRC=[];
M_values=[1,2,4,8];
for m=1:6% for each experiment
    for target=0.1:0.1:10
    Pout_SC=0;
    Pout_MRC=0;
    
        for t=1:500
        M=M_values(m);% # of branches (independent fading paths)
        max_snr_SC=0;
        max_snr_MRC=0;
        
        for j=1:M
            
            %complex signal with gaussian noise
            sig = 1i*ones(2000,1)*(j/10)+10^-3*wgn(2000,1,0);        
            out = rayChan(sig);
            %plot(20*log10(abs(out)))
            
            alpha=10/j;%coefficients
            out=alpha*abs(out);%cophasing
            beta=snr(out);%take the SNR
            
            %selection, SC
            if(max_snr_SC<beta)
                max_snr_SC=beta;
            end
            
            %selection, MRC
            %sum all SNRs not focus on just one branch
            max_snr_MRC=max_snr_MRC+beta;
        end
        SC_snr=10^(max_snr_SC/20);%convert db to watt
        MRC_snr=10^(max_snr_MRC/20);%convert db to watt
        
        %update Pout
        if SC_snr<target
           Pout_SC=Pout_SC+1;
        end
        
        if MRC_snr<target
           Pout_MRC=Pout_MRC+1;
        end
        
        end
    outage_values_SC(m,round(10*target))=Pout_SC/t;
    outage_values_MRC(m,round(10*target))=Pout_MRC/t;
    end
   
end