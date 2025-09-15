%   Calculating the SER and BER of M-QAM

M = 8;
k = log2(M);
K = 10; % Rician K ratio index
EbNodB = -30:30;
h = getRician(K,length(EbNodB));
EbNodBFaded = 10*log10(abs(h).^2.*10.^((EbNodB)./10));
snrdB = convertSNR(EbNodBFaded,"ebno","snr",BitsPerSymbol=k);
% instSNR = abs(h).^2.*snrdB;
g = mean(EbNodB);
figure;
hold on;
for i=1:length(M)
    berTheoretical = berfading(EbNodB,"qam",M(i),1,10);
    a = 0.2*(1+K)*(M(i)-1)*log2(M(i));
    b = (1+K)*(M(i)-1)+1.5*10.^(EbNodBFaded./10);
    c = 1.5*10.^(EbNodBFaded./10)*K./b;
    Ps = (a./b).*exp(-c);
    semilogy(EbNodB,berTheoretical,'-*');
    semilogy(EbNodBFaded,Ps./k(i),'-o');
end
hold off;
title('BER Comparison in Rician Fading for Different Modulation Orders');
xlabel('Eb/No [dB]');
ylabel('Bit-Error Rate (BER)');
legend('berfading','Paper by Ali et. al');
grid on;

function h = getRician(K,N)
%   Draw one rician coefficient from the distribution
    hLOS = exp(1i*2*pi*rand(1,N));
    hNLOS = sqrt(1/2)*(randn(1, N) + 1i*randn(1, N));
    h = sqrt(K/(K + 1))*hLOS + sqrt(1/(K + 1))*hNLOS;
end