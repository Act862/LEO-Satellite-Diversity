%% 4x1 MISO with Rician fading, delays, guard spacing, and MRC combining
% Fixes overlapping-delays issue by upsampling (inserting zeros between symbols)
clear; clc; rng(2025);

%% Parameters
msg        = 'Hello';                    % message to send
bits       = reshape(de2bi(double(msg),8,'left-msb').',[],1); % column vector of bits
Nbits      = length(bits);               % 5 chars * 8 = 40 bits
bpsk       = 2*bits - 1;                 % BPSK: 0->-1, 1->+1 (real)
Nt         = 4;                          % # transmit antennas
delays     = [0 2 4 8];                  % delays in *samples* (integer)
maxDelay   = max(delays);
SNRdB      = 20;                         % SNR in dB (adjustable)
SNRlin     = 10^(SNRdB/10);

%% Upsampling / guard spacing
% Choose factor L so that symbols are spaced > maxDelay to avoid overlap
L = maxDelay + 1;         % simplest safe choice: spacing is at least maxDelay
% Map each symbol to one filled sample followed by (L-1) zeros
Ns_symbols = Nbits;
Ns = Ns_symbols * L;      % total number of samples transmitted from each antenna

tx_base = zeros(Ns,1);
tx_base(1:L:Ns) = bpsk;   % place symbol every L samples

% Build per-antenna transmit signals by shifting according to delays
txSig = zeros(Ns + maxDelay, Nt); % extra tail to hold later shifted samples
for t = 1:Nt
    shift = delays(t);
    txSig(1+shift : 1+shift + Ns -1, t) = tx_base;
end
Ns_total = size(txSig,1);

%% Channel: Rician per-antenna (flat per-sample), same across time for simplicity
K = 5; % Rician K-factor
h_LOS = ones(1,Nt); % deterministic LOS phasors (choose 1 for simplicity)
h_NLOS = (randn(1,Nt)+1j*randn(1,Nt))/sqrt(2); % diffuse
h = sqrt(K/(K+1))*h_LOS + sqrt(1/(K+1))*h_NLOS; % 1 x Nt complex gains

% Apply channel (linear): y[n] = sum_t h(t) * x_t[n] + noise[n]
y = zeros(Ns_total,1);
for t = 1:Nt
    y = y + h(t) * txSig(:,t);
end

% Add AWGN (complex) with energy normalization
Es_symbol = mean(abs(tx_base(tx_base~=0)).^2); % should be 1
N0 = Es_symbol / SNRlin;
noise = sqrt(N0/2)*(randn(Ns_total,1) + 1j*randn(Ns_total,1));
r = y + noise;

%% Receiver: time-align each branch, apply MRC (conjugate weighting) and sample
% Extract per-branch contributions by shifting the received sequence by -delay
% For each t: the portion of r corresponding to tx_base is located at indices
%   idx = (1:Ns) + delays(t)
% So to align to symbol positions (1:L:Ns), we extract r(1+delays(t) : L : 1+delays(t)+Ns-1)
symbol_positions = 1:L:Ns; % indices in tx_base containing symbol pulses (length Nbits)
rx_per_branch = zeros(Nbits, Nt);

for t = 1:Nt
    idx_start = 1 + delays(t);         % first sample index in r for this antenna's first symbol
    % indices in r that carry the symbol samples for this antenna
    idxs = idx_start : L : idx_start + (Nbits-1)*L;
    % ensure indices are within bounds
    if idxs(end) > length(r)
        error('Indexing exceeded received vector length - check padding/lengths.');
    end
    rx_samples_t = r(idxs).';            % 1 x Nbits complex (transpose)
    % MRC weight = conj(h(t)) to undo phase of channel
    rx_per_branch(:, t) = (conj(h(t)) * rx_samples_t).'; % Nbits x 1
end

% Combine branches (sum across antennas)
r_combined = sum(rx_per_branch, 2); % Nbits x 1 complex

% Decision for BPSK: symbols are real +/-1 after MRC (phase removed)
bits_hat = real(r_combined) > 0;  % logical vector length Nbits

% Reconstruct message
msg_hat = char(bi2de(reshape(bits_hat, 8, []).', 'left-msb')).';

%% Results
disp(['Original Message: ' msg]);
disp(['Decoded Message : ' msg_hat]);
numErr = sum(bits ~= bits_hat);
disp(['Bit errors: ' num2str(numErr) ' / ' num2str(Nbits)]);
disp('Channel gains (h):');
disp(h);

% (Optional) show a small constellation snapshot (real vs imag of combined samples)
figure;
plot(real(r_combined), imag(r_combined), 'o');
xlabel('Real'); ylabel('Imag'); grid on;
title('Received symbols after MRC (each sample corresponds to a BPSK symbol)');
