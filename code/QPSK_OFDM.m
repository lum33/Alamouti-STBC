% QPSK_OFDM.m
% Compare 1x1, 2x1, 2x2 Alamouti STBC (QPSK) + OFDM
%% Parameter
L_frame = 128;      % Number of Alamouti blocks/packet
N_Packets = 10000;  % Number of packets
NT = 2; 
NR = 1:1:2;                 
M = 4;              % QPSK
b = log2(M); 
offset = pi/4;      % Phase offset for QPSK symbols  
SNRdBs = -5:2:25;

% OFDM Parameters
N_subcarriers = 64; 
CP_length = 32;    
N_OFDM = L_frame / N_subcarriers;

%% 1. QPSK Simulation (no STBC + no OFDM)
figure(1)
BER_qpsk = zeros(length(SNRdBs), 1);
total_bits_qpsk = L_frame * b;

for idx = 1:length(SNRdBs)
    SNRdB = SNRdBs(idx);
    sigma = sqrt(1/(2*10^(SNRdB/10))); % noise std
    
    err_count = 0;
    
    for packet = 1:N_Packets
        % Generate random data
        bits = randi([0 1], total_bits_qpsk, 1);
        X_int = bi2de(reshape(bits, b, []).', 'left-msb');
        % Modulation: QPSK
        X_mod = pskmod(X_int, M, offset);
        
        H = (randn(1,1) + 1j*randn(1,1)) / sqrt(2); % (1 x 1)
        Y = H * X_mod + sigma*(randn(1, L_frame) + 1j*randn(1, L_frame));
        Y_eq = Y / H; % Equalization
        
        % Maximum Likelihood Detection
        refSymbols = pskmod((0:M-1).', M, offset);
        s_hat_ml = zeros(L_frame,1);
        for n = 1:L_frame
            y = Y_eq(n);
            distances = abs(y - refSymbols).^2;
            [~, minIdx] = min(distances);
            s_hat_ml(n) = refSymbols(minIdx);
        end
        
        % Demodulation
        X_hat = pskdemod(s_hat_ml, M, offset);
        err_count = err_count + sum(X_hat ~= X_int);
    end
    BER_qpsk(idx) = err_count / (total_bits_qpsk * N_Packets);
end

semilogy(SNRdBs, BER_qpsk, 'g-^', 'LineWidth', 2);  % QPSK curve
grid on;
hold on;

%% 2. QPSK Simulation (no STBC + OFDM)
BER_qpsk_OFDM = zeros(length(SNRdBs), 1);

for idx = 1:length(SNRdBs)
    SNRdB = SNRdBs(idx);
    sigma = sqrt(1/(2*10^(SNRdB/10))); % Noise standard deviation
    
    err_count = 0;
    
    for packet = 1:N_Packets
        % Generate random data
        bits = randi([0 1], total_bits_qpsk, 1);
        X_int = bi2de(reshape(bits, b, []).', 'left-msb');

        % Modulation: QPSK
        X_mod = pskmod(X_int, M, offset);
        X_mod_matrix = reshape(X_mod, N_subcarriers, N_OFDM);

        tx_OFDM = zeros(N_subcarriers + CP_length, N_OFDM);
        for k = 1:N_OFDM
            % Apply IFFT
            ifft_data = ifft(X_mod_matrix(:,k), N_subcarriers);
            % Add CP
            tx_OFDM(:,k) = [ifft_data(end-CP_length+1:end); ifft_data];
        end

        tx_signal = tx_OFDM(:).';
        tx_signal = tx_signal / sqrt(mean(abs(tx_signal).^2));
        
        H = (randn(1,1) + 1j*randn(1,1)) / sqrt(2); % (1 x 1)
        rx_signal = H * tx_signal + sigma * (randn(size(tx_signal)) + 1j*randn(size(tx_signal)));
        
        rx_OFDM = reshape(rx_signal, N_subcarriers + CP_length, N_OFDM);
        Y = zeros(N_subcarriers, N_OFDM);
        for k = 1:N_OFDM
            % Remove CP
            rx_noCP = rx_OFDM(CP_length+1:end, k);
            % Apply FFT
            Y(:, k) = fft(rx_noCP, N_subcarriers);
        end

        % Equalization
        Y_eq = Y / H;
        Y_eq_vec = Y_eq(:);

        % Maximum Likelihood Detection
        refSymbols = pskmod((0:M-1).', M, offset);
        s_hat_ml = zeros(length(Y_eq_vec), 1);
        for n = 1:length(Y_eq_vec)
            y = Y_eq_vec(n);
            distances = abs(y - refSymbols).^2;
            [~, minIdx] = min(distances);
            s_hat_ml(n) = refSymbols(minIdx);
        end
        
        % Demodulation
        X_hat = pskdemod(s_hat_ml, M, offset);
        err_count = err_count + sum(X_hat ~= X_int);
    end
    BER_qpsk_OFDM(idx) = err_count / (total_bits_qpsk * N_Packets);
end

% Plot BER
semilogy(SNRdBs, BER_qpsk_OFDM, 'b-o', 'LineWidth', 2); % OFDM-QPSK curve
grid on;
hold on;
legend('SISO (1Tx, 1Rx)', 'OFDM-SISO (1Tx, 1Rx)')
l=legend(BackgroundAlpha=.5);
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
title('BER vs SNR for QPSK with OFDM');
saveas(figure(1),'BER_OFDM_comparison.jpg')