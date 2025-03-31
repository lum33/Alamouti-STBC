%% 3. 2x1, 2x2 Alamouti STBC (no OFDM) 
BER_STBC = zeros(length(SNRdBs), length(NR));

markers = {'o', 's'};
colors = {'r', 'b'};

% Switch between 2x1 and 2x2 case
for r = 1:length(NR)
    % Simulation over SNR values
    for idx = 1:length(SNRdBs)
        SNRdB = SNRdBs(idx);
        sigma = sqrt(1/(2*10^(SNRdB/10))); % noise std
    
        err_count = 0;
    
        for packet = 1:N_Packets
            % Generate Random Data
            bits = randi([0 1], L_frame * b, 1);
            X_int = bi2de(reshape(bits, b, []).', 'left-msb');
            X_mod = pskmod(X_int, M, offset);
            
            % Alamouti Encoding
            [tx_1, tx_2] = AlamoutiEncoder(X_mod, NT);
            
            % Channel Simulation
            H = (randn(r, NT) + 1j*randn(r, NT)) / sqrt(NT); % (NR x NT)
            
            % H for decoding (for each symbol block)
            H_block = repmat(H, [1, 1, L_frame]);
            
            % Transmit through the channel over two time slots
            rx1_signal = zeros(r, L_frame);
            rx2_signal = zeros(r, L_frame);
            for rx = 1:r
                rx1_signal(rx,:) = sum(H(rx,:) * tx_1, 1) + sigma*(randn(1,L_frame) + 1j*randn(1,L_frame));
                rx2_signal(rx,:) = sum(H(rx,:) * tx_2, 1) + sigma*(randn(1,L_frame) + 1j*randn(1,L_frame));
            end
            
            % Alamouti Decoding
            s_hat = AlamoutiDecoder(rx1_signal, rx2_signal, H_block);

            % Maximum Likelihood Detection
            refSymbols = pskmod((0:M-1).', M, offset); % reference QPSK constellation points

            s_hat_ml = zeros(size(s_hat));
            for tx = 1:NT
                for n = 1:L_frame
                    y = s_hat(tx, n);
                    distances = abs(y - refSymbols).^2;  % Euclidean distance squared
                    [~, minIdx] = min(distances);
                    s_hat_ml(tx, n) = refSymbols(minIdx);
                end
            end
            
            % Demodulation
            X_hat_int_ml = pskdemod(s_hat_ml, M, offset);
            
            % Count symbol errors for each transmit antenna
            for tx = 1:NT
                err_count = err_count + sum(X_hat_int_ml(tx,:) ~= X_int(tx,:));
            end
        end
        % BER calculation
        BER_STBC(idx, r) = err_count / (total_bits * N_Packets);
    end
    % Plot BER vs SNR
    semilogy(SNRdBs, BER_STBC(:, r), ...
            [colors{r} '-' markers{r}], ...
            'LineWidth', 2);
    grid on;
    hold on;

end

%% 4. 2x1, 2x2 Alamouti STBC (OFDM) 
BER_STBC_OFDM = zeros(length(SNRdBs), length(NR));

markers_OFDM = {'o', 's'};
colors_OFDM = {'r', 'b'};

% Switch between 2x1 and 2x2 case
for r = 1:length(NR)
    % Simulation over SNR values
    for idx = 1:length(SNRdBs)
        SNRdB = SNRdBs(idx);
        sigma = sqrt(1/(2*10^(SNRdB/10))); % noise std
    
        err_count = 0;
    
        for packet = 1:N_Packets
            % Generate Random Data
            bits = randi([0 1], L_frame * b, 1);
            X_int = bi2de(reshape(bits, b, []).', 'left-msb');
            X_mod = pskmod(X_int, M, offset);
            
            % Alamouti Encoding
            [tx_1, tx_2] = AlamoutiEncoder(X_mod, NT);

            tx_1_matrix = reshape(tx_1, N_subcarriers, N_OFDM);
            tx_2_matrix = reshape(tx_2, N_subcarriers, N_OFDM);
            
            tx1_OFDM = zeros(N_subcarriers + CP_length, N_OFDM);
            tx2_OFDM = zeros(N_subcarriers + CP_length, N_OFDM);
            for k = 1:N_OFDM
                % Apply IFFT
                ifft_data_1 = ifft(tx_1_matrix(:,k), N_subcarriers);
                ifft_data_2 = ifft(tx_2_matrix(:,k), N_subcarriers);
                % Add CP
                tx1_OFDM(:,k) = [ifft_data_1(end-CP_length+1:end); ifft_data_1];
                tx2_OFDM(:,k) = [ifft_data_2(end-CP_length+1:end); ifft_data_2];
            end
    
            tx1_signal = tx1_OFDM(:).';
            tx1_signal = tx1_signal / sqrt(mean(abs(tx1_signal).^2));
            tx2_signal = tx2_OFDM(:).';
            tx2_signal = tx2_signal / sqrt(mean(abs(tx2_signal).^2));
            
            % Channel Simulation
            H = (randn(r, NT) + 1j*randn(r, NT)) / sqrt(NT); % (NR x NT)
            
            % H for decoding (for each symbol block)
            H_block = repmat(H, [1, 1, L_frame]);
            
            % Transmit through the channel over two time slots
            rx1_signal = zeros(r, L_frame);
            rx2_signal = zeros(r, L_frame);
            for rx = 1:r
                rx1_signal(rx,:) = sum(H(rx,:) * tx_1, 1) + sigma*(randn(1,L_frame) + 1j*randn(1,L_frame));
                rx2_signal(rx,:) = sum(H(rx,:) * tx_2, 1) + sigma*(randn(1,L_frame) + 1j*randn(1,L_frame));
            end

            rx1_OFDM = reshape(rx1_signal, N_subcarriers + CP_length, N_OFDM);
            rx2_OFDM = reshape(rx2_signal, N_subcarriers + CP_length, N_OFDM);
            Y1 = zeros(N_subcarriers, N_OFDM);
            Y2 = zeros(N_subcarriers, N_OFDM);
            for k = 1:N_OFDM
                % Remove CP
                rx1_noCP = rx1_OFDM(CP_length+1:end, k);
                rx2_noCP = rx2_OFDM(CP_length+1:end, k);
                % Apply FFT
                Y1(:, k) = fft(rx1_noCP, N_subcarriers);
                Y2(:, k) = fft(rx2_noCP, N_subcarriers);
            end

            % Equalization
            Y1_eq = Y1 / H;
            Y2_eq = Y2 / H;
            Y1_eq_vec = Y1_eq(:);
            Y2_eq_vec = Y2_eq(:);


            % Alamouti Decoding
            s_hat = AlamoutiDecoder(Y1_eq_vec, Y2_eq_vec, H_block);

            % Maximum Likelihood Detection
            refSymbols = pskmod((0:M-1).', M, offset); % reference QPSK constellation points

            s_hat_ml = zeros(size(s_hat));
            for tx = 1:NT
                for n = 1:L_frame
                    y = s_hat(tx, n);
                    distances = abs(y - refSymbols).^2;  % Euclidean distance squared
                    [~, minIdx] = min(distances);
                    s_hat_ml(tx, n) = refSymbols(minIdx);
                end
            end
            
            % Demodulation
            X_hat_int_ml = pskdemod(s_hat_ml, M, offset);
            
            % Count symbol errors for each transmit antenna
            for tx = 1:NT
                err_count = err_count + sum(X_hat_int_ml(tx,:) ~= X_int(tx,:));
            end
        end
        % BER calculation
        BER_STBC_OFDM(idx, r) = err_count / (total_bits * N_Packets);
    end
    % Plot BER vs SNR
    semilogy(SNRdBs, BER_STBC_OFDM(:, r), ...
            [colors_OFDM{r} '-' markers_OFDM{r}], ...
            'LineWidth', 2);
    grid on;
    hold on;

end

legend('Alamouti (2Tx, 1Rx)', 'Alamouti (2Tx, 2Rx)',...
    'Alamouti-OFDM (2Tx, 1Rx)', 'Alamouti-OFDM (2Tx, 2Rx)')
l=legend(BackgroundAlpha=.5);
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
title('BER vs SNR for QPSK with OFDM');
saveas(figure(1),'BER_STBC_OFDM.jpg')