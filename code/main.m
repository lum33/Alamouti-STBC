% main.m
%% Parameters
L_frame = 130;                  % Number of Alamouti blocks/packet (each block carries 2 symbols)
N_Packets = 4000;               % Number of packets
NT = 2; 
NR = 2;                 
M = 4;                          % QPSK
b = log2(M); 
total_bits = NT * L_frame * b;
offset = pi/4;                  % Phase offset for QPSK symbols  
SNRdBs = 0:2:30;
BER = zeros(length(SNRdBs), 1);

% GIF parameter
gifFileName = 'RxConstellation_2X2.gif';
gifFrameDelay = 0.01; % seconds delay between frames
hGifFig = figure('Name', 'Receiver Constellation for GIF');

%% Simulation over SNR values
for idx = 1:length(SNRdBs)
    SNRdB = SNRdBs(idx);
    sigma = sqrt(1/(2*10^(SNRdB/10))); % noise std

    err_count = 0;

    for packet = 1:N_Packets
        %% Generate Random Data and QPSK Modulation
        bits = randi([0 1], total_bits, 1);
        % Map bits to symbols for each tx antenna
        X_int = zeros(NT, L_frame);  % integer symbols [0, M-1]
        for tx = 1:NT
            bits_tx = bits(((tx-1)*L_frame*b + 1):(tx*L_frame*b)); % Each antenna gets its portion of bits
            X_int(tx,:) = bi2de(reshape(bits_tx, b, []).', 'left-msb');
        end
        
        % Modulation: QPSK
        X_mod = pskmod(X_int, M, offset);
        
        %% Alamouti Encoding
        [tx_1, tx_2] = AlamoutiEncoder(X_mod);
        
        %% Channel Simulation
        H = (randn(NR, NT) + 1j*randn(NR, NT)) / sqrt(2); % (NR x NT)
        
        % H for decoding (for each symbol block)
        H_block = repmat(H, [1, 1, L_frame]);
        
        % Transmit through the channel over two time slots
        Y1 = zeros(NR, L_frame);
        Y2 = zeros(NR, L_frame);
        for rx = 1:NR
            Y1(rx,:) = sum(H(rx,:) * tx_1, 1) + sigma*(randn(1,L_frame) + 1j*randn(1,L_frame));
            Y2(rx,:) = sum(H(rx,:) * tx_2, 1) + sigma*(randn(1,L_frame) + 1j*randn(1,L_frame));
        end
        
        %% Alamouti Decoding
        s_hat = AlamoutiDecoder(Y1, Y2, H_block);
        
        %% Maximum Likelihood Detection
        refSymbols = pskmod((0:M-1).', M, offset); % reference QPSK constellation points
        
        % ML detection
        s_hat_ml = zeros(size(s_hat));
        for tx = 1:NT
            for n = 1:L_frame
                y = s_hat(tx, n);
                distances = abs(y - refSymbols).^2;  % Euclidean distance squared
                [~, minIdx] = min(distances);
                s_hat_ml(tx, n) = refSymbols(minIdx);
            end
        end
        
        %% Demodulation and BER Calculation
        X_hat_int_ml = pskdemod(s_hat_ml, M, offset);
        
        % Count symbol errors for each transmit antenna
        for tx = 1:NT
            err_count = err_count + sum(X_hat_int_ml(tx,:) ~= X_int(tx,:));
        end

        %% Record Receiver Constellation for GIF 
        % Concatenate the decoded symbols for visualization
        if packet >= 1 && packet <= 20
            rx_symbols = [s_hat(1,:); s_hat(2,:)];
            figure(hGifFig);
            scatter(real(rx_symbols), imag(rx_symbols), 36, 'r', 'filled');
            axis([-1.5 1.5 -1.5 1.5]);
            grid on;
            title(sprintf('Receiver Constellation (Packet %d, SNR = %d dB)', packet, SNRdB));
            drawnow;
            
            % Capture the current figure as a frame and write to GIF
            frame = getframe(hGifFig);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im,256);
            if idx == 1
                imwrite(imind, cm, gifFileName, 'gif', 'Loopcount', inf, 'DelayTime', gifFrameDelay);
            else
                imwrite(imind, cm, gifFileName, 'gif', 'WriteMode', 'append', 'DelayTime', gifFrameDelay);
            end
        end
    end
    % BER calculation
    BER(idx) = err_count / (total_bits * N_Packets);
end

%% Plot BER vs SNR
figure(1);
semilogy(SNRdBs, BER, 'b-o', 'LineWidth', 2);
grid on;
xlabel('SNR [dB]');
ylabel('Bit Error Rate (BER)');
%axis([SNRdBs([1 end]) 1e-6 1e0]); 
title('BER vs SNR for 2x2 Alamouti Scheme with QPSK');
saveas(figure(1),'BER_2X2.jpg')