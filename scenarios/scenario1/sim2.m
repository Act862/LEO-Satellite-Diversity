%   Implementing several diversity combining strategies
%   L = (1,2,3,4)
%   SIMO - Single Input Multiple Output
%   Uncorrelated Rayleigh fading
%   Average symbol energy-to-noise power ratio Es/N0 = 9dB
%   Simulate QPSK error rate for
%   1. Selective Combining
%   2. Maximal Ratio Combining
%   3. Equal Gain Combining
%   4. Direct Combining
sample_num = 1e6;
enr_range = 0:16;
Pe_mrc = zeros(length(enr_range),4);
%   Rayleigh Fading
sigma = 1/sqrt(2);
for enr_index = 1:length(enr_range)
    enr_dB = enr_range(enr_index); % energy to noise ratio
    enr = 10^(enr_dB/10);
    %   generate QPSK data (2 bits)
    data = rand(2,sample_num);
    data = 2*(data > 0.5)-1;    % map to -1,1
    Edata = sqrt(2); % symbol energy
    En = Edata/enr; % noise energy
    %   for each channel
    for L=1:4
        % n = normrnd(0,sqrt(En/2),2,sample_num,L) + ...
        %     1i*normrnd(0,sqrt(En/2),2,sample_num,L);
        n = randn(2,sample_num,L)*sqrt(En/2) + 1i*randn(2,sample_num,L);
        %   generate fading gain
        g = (randn(1,sample_num,L)*sigma) + 1i*(1/2*randn(1,sample_num,L));
        %   repeat the gain for each channel
        g_tmp = repmat(g,2,1,1);
        tx_data = repmat(data,1,1,L);
        r = g_tmp.*tx_data + n;
        % Maximal Ratio Combining
        [Pe_mrc(enr_index,L), result_mrc] = mrc(g_tmp,r,sample_num,data);
    end
end

figure(1),semilogy(enr_range, Pe_mrc(:,:),'-*')
grid on;
title('BER of Maximal Ratio Combining (Rayleigh)');
legend('L=1','L=2','L=3','L=4');
xlabel('SNR (dB)');
ylabel('Bit error rate');
