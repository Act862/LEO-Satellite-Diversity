% Parameters
num_trials = 1000;
M_values = [1, 4]; % diversity levels
snr_threshold_dB = 5;
P_tx_dBm = -80; % transmit power
N0_dBm = -100; % noise floor
snr_linear = 10.^((P_tx_dBm - N0_dBm)/10);
snr_threshold = 10^(snr_threshold_dB/10);

% Delay settings
base_delay = 15e-3; % 15 ms fixed base delay from earlier sim
retransmit_penalty = 5e-3; % 5 ms per retransmission

% Storage
avg_delay = zeros(size(M_values));
packet_loss_prob = zeros(size(M_values));

for m = 1:length(M_values)
    M = M_values(m);
    total_delay = 0;
    failed_packets = 0;

    for n = 1:num_trials
        attempts = 0;
        success = false;

        while attempts < 4 && ~success
            % Generate log-logistic fading samples
            fading_samples = lognrnd(0,1,1,M); % log-logistic ~ log-normal approx.
            % Diversity combining (Selection combining = max)
            combined_snr = max(fading_samples) * snr_linear;

            if combined_snr >= snr_threshold
                success = true;
                total_delay = total_delay + base_delay + retransmit_penalty * attempts;
            else
                attempts = attempts + 1;
            end
        end

        if ~success
            failed_packets = failed_packets + 1;
            total_delay = total_delay + base_delay + retransmit_penalty * 3; % max penalty
        end
    end

    avg_delay(m) = total_delay / num_trials;
    packet_loss_prob(m) = failed_packets / num_trials;
end

% Plotting
figure;
bar(avg_delay * 1e3);
set(gca, 'XTickLabel', {'M=1', 'M=4'});
ylabel('Average Delay (ms)');
title('Effect of Spatial Diversity on Delay with Retransmissions');
grid on;

figure;
bar(packet_loss_prob * 100);
set(gca, 'XTickLabel', {'M=1', 'M=4'});
ylabel('Packet Loss Probability (%)');
title('Effect of Spatial Diversity on Packet Loss (Log-Logistic Fading)');
grid on;
