x = 1; % signal to transmit Eb = 1
TRIAL = 10000; %number of simulation runs per EbN0 %50000
for EbN0 = 0:1:20 %dB
    linear_EbN0 = 10^(EbN0/10); nvar = 1/(linear_EbN0); %calculation of N0, remember Eb = 1
    error1 = 0; %set error counter to 0
    error2 = 0; %set error counter to 0
    error3 = 0; %set error counter to 0
        for trial = 1:TRIAL % monte carlo trials.. count the errors
            n1 = sqrt(nvar/2)*randn; %noise for the first
            n2 = sqrt(nvar/2)*randn; %noise for the first
            h1 = sqrt(0.5)*abs(randn + j*randn); %rayleigh amplitude 1
            h2 = sqrt(0.5)*abs(randn + j*randn); %rayleigh amplitude 1

            %Equal Gain combining
            y1 = x*h1+n1; % Signal 1
            y2 = x*h2+n2; % Signal 2
            y_equal = 0.5*(y1+y2); 

            %Maximal Ratio combining
            a1 = (abs(h1))^2;
            a2 = (abs(h2))^2;
            y_maximal = x*(a1*h1+a2*h2)+a1*n1+a2*n2;

            %Selection combining
            P1 = chi2rnd(4);
            P2 = chi2rnd(4);
            as1 = P1*(abs(h1))^2;
            as2 = P2*(abs(h2))^2;
            if as1 >= as2
                y_selection = x*(as1*h1)+as1*n1;
            end
            if as1 < as2
                y_selection = x*(as2*h2)+as2*n2;
            end

            if y_equal < 0 %define decision region as 0 
                error1 = error1 + 1;
            end
            if y_maximal < 0 
                error2 = error2 + 1;
            end
            if y_selection < 0 
                error3 = error2 + 1;
            end
        end
    BER1(EbN0+1) = error1/(TRIAL);
    BER2(EbN0+1) = error2/(TRIAL);
    BER3(EbN0+1) = error3/(TRIAL);
end
% plot simulations
figure
EbNo=0:1:20; %changed from 10
mu = 10.^(EbNo./10);
ber_theory = (1/2)*(1 - sqrt(mu ./ (mu + 1))); 
semilogy(EbNo,BER1,'r*-',EbNo,BER2,'b--o',EbNo,BER3,'c-o',EbNo,ber_theory,'b'); % plot EG BER vs EbNo 
legend('EG','MR','SC','theory');
xlabel('EbNo(dB)') %Label for x-axis
ylabel('Bit error rate') %Label for y-axis