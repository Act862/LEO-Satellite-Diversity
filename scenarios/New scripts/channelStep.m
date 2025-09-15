function [nextState, reward] = channelStep()
    % --- simulate channel as before ---
    K = 5;  % Rician K-factor
    phi = 2*pi*rand;
    LOS = sqrt(K/(K+1));
    sigma = sqrt(1/(2*(K+1)));
    h_complex = LOS*exp(1j*phi) + sigma*(randn + 1j*randn);
    h_abs2 = abs(h_complex)^2;

    % channel estimation
    pilotSNR_dB = 10; Npilots = 4;
    pilotSNR_lin = 10^(pilotSNR_dB/10);
    sigma_e2 = 1/(pilotSNR_lin * Npilots);
    e_complex = sqrt(sigma_e2/2)*(randn + 1j*randn);
    h_hat = h_complex + e_complex;
    U = sigma_e2 / (abs(h_hat)^2 + eps);

    % nominal Eb/N0
    ebno_nom_dB = 10*rand;  % between 0 and 10 dB
    ebno_nom_lin = 10^(ebno_nom_dB/10);

    % compute true SNRs
    ebno_faded_lin = ebno_nom_lin * h_abs2;
    snr_sc = 10*log10(ebno_faded_lin);

    ebno_est_faded_lin = ebno_nom_lin * abs(h_hat)^2;
    snr_mrc = 10*log10(ebno_est_faded_lin);

    % package state
    nextState = [U, abs(h_hat)^2, ebno_nom_dB];

    % store both possible rewards so agent can pick later
    reward = [snr_sc, snr_mrc];
end
