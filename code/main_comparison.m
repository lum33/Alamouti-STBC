% main_comparison.m
% Compare 1x1, 2x1, 2x2 Alamouti STBC (QPSK) 
%% 2x1, 2x2 Alamouti STBC (QPSK) 
L_frame = 130;                  % Number of Alamouti blocks/packet (each block carries 2 symbols)
N_Packets = 4000;               % Number of packets
NT = 2; 
NR = 1:1:2;                 
M = 4;                          % QPSK
b = log2(M); 
total_bits = NT * L_frame * b;
offset = pi/4;                  % Phase offset for QPSK symbols  
SNRdBs = -15:2:20;
BER = zeros(length(SNRdBs), length(NR));

figure(1)
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
            bits = randi([0 1], total_bits, 1);
            % Map bits to symbols for each tx antenna
            X_int = zeros(NT, L_frame);  % integer symbols [0, M-1]
            X_int_QPSK = zeros(1, L_frame);

            for tx = 1:NT
                bits_tx = bits(((tx-1)*L_frame*b + 1):(tx*L_frame*b)); % Each antenna gets its portion of bits
                X_int(tx,:) = bi2de(reshape(bits_tx, b, []).', 'left-msb');
            end
            
            % Modulation: QPSK
            X_mod = pskmod(X_int, M, offset);
            
            % Alamouti Encoding
            [tx_1, tx_2] = AlamoutiEncoder(X_mod);
            
            % Channel Simulation
            H = (randn(r, NT) + 1j*randn(r, NT)) / sqrt(2); % (NR x NT)
            
            % H for decoding (for each symbol block)
            H_block = repmat(H, [1, 1, L_frame]);
            
            % Transmit through the channel over two time slots
            Y1 = zeros(r, L_frame);
            Y2 = zeros(r, L_frame);
            for rx = 1:r
                Y1(rx,:) = sum(H(rx,:) * tx_1, 1) + sigma*(randn(1,L_frame) + 1j*randn(1,L_frame));
                Y2(rx,:) = sum(H(rx,:) * tx_2, 1) + sigma*(randn(1,L_frame) + 1j*randn(1,L_frame));
            end
            
            % Alamouti Decoding
            s_hat = AlamoutiDecoder(Y1, Y2, H_block);
            
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
        BER(idx, r) = err_count / (total_bits * N_Packets);
    end
    % Plot BER vs SNR
    semilogy(SNRdBs, BER(:, r), ...
            [colors{r} '-' markers{r}], ...
            'LineWidth', 2);
    grid on;
    hold on;

end

%% QPSK Simulation (no STBC)
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

semilogy(SNRdBs, BER_qpsk, 'g-^', 'LineWidth', 2);  % QPSK curve in black dashed line with triangle markers

legend('Alamouti (2x1)', 'Alamouti (2x2)', 'QPSK (1x1)')
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
%axis([SNRdBs([1 end]) 1e-6 1e0]);
title('BER vs SNR with QPSK');
saveas(figure(1),'BER_comparison.jpg')