%% Alamouti_2x2_QPSK_GIF.m
clear; close all; clc;

%% Parameters
L_frame = 130;          % Number of Alamouti blocks/packet (each block carries 2 symbols)
N_Packets = 4000;       % Number of packets
NT = 2; 
NR = 2;                 
M = 4;                  % QPSK
b = log2(M);            
SNRdBs = 0:2:30;        
offset = pi/4;      % Phase offset so that the constellation lies on the four corners of a square
txConstDiagram = comm.ConstellationDiagram('Title', 'Transmitter Constellation Diagram');
BER = zeros(length(SNRdBs), 1);

% GIF parameter
gifFileName = 'rx_constellation.gif';
gifFrameDelay = 0.01; % seconds delay between frames
hGifFig = figure('Name', 'Receiver Constellation for GIF');

%% Simulation over SNR values
for idx = 1:length(SNRdBs)
    SNRdB = SNRdBs(idx);
    sigma = sqrt(1/(2*10^(SNRdB/10))); % noise std
    
    err_count = 0;
    total_bits = 0;
    
    for packet = 1:N_Packets
        %% Transmitter: Generate and Modulate Data
        % Generate random bits for both streams
        bits = randi([0 1], 2*L_frame*b, 1);
        % Map bits to symbols (each symbol uses b bits)
        data_symbols = bi2de(reshape(bits, b, []).', 'left-msb');
        
        % Split data into two streams
        s1 = data_symbols(1:L_frame);
        s2 = data_symbols(L_frame+1:end);
        
        % Modulation: QPSK
        x1 = pskmod(s1, M, offset);
        x2 = pskmod(s2, M, offset);
        
        % 2x2 independent Rayleigh Fading Channel        
        h1_1 = (randn(L_frame,1) + 1j*randn(L_frame,1))/sqrt(2);
        h2_1 = (randn(L_frame,1) + 1j*randn(L_frame,1))/sqrt(2);
        h1_2 = (randn(L_frame,1) + 1j*randn(L_frame,1))/sqrt(2);
        h2_2 = (randn(L_frame,1) + 1j*randn(L_frame,1))/sqrt(2);
        
        % Received signals with Alamouti encoding
        r1_1 = h1_1.*x1 + h2_1.*x2;
        r2_1 = h1_1.*(-conj(x2)) + h2_1.*(conj(x1));
        r1_2 = h1_2.*x1 + h2_2.*x2;
        r2_2 = h1_2.*(-conj(x2)) + h2_2.*(conj(x1));
        
        %% Receiver
        r1_1 = r1_1 + sigma*(randn(L_frame,1) + 1j*randn(L_frame,1));
        r2_1 = r2_1 + sigma*(randn(L_frame,1) + 1j*randn(L_frame,1));
        r1_2 = r1_2 + sigma*(randn(L_frame,1) + 1j*randn(L_frame,1));
        r2_2 = r2_2 + sigma*(randn(L_frame,1) + 1j*randn(L_frame,1));
        
        % Alamouti Decoding
        y1 = (r1_1.*conj(h1_1) + conj(r2_1).*h2_1) + (r1_2.*conj(h1_2) + conj(r2_2).*h2_2);
        y2 = (r1_1.*conj(h2_1) - conj(r2_1).*h1_1) + (r1_2.*conj(h2_2) - conj(r2_2).*h1_2);
        
        % Normalize with total channel power
        h_power = (abs(h1_1).^2 + abs(h2_1).^2) + (abs(h1_2).^2 + abs(h2_2).^2); 
        y1 = y1 ./ h_power;
        y2 = y2 ./ h_power;
        
        % Demodulation
        s1_hat = pskdemod(y1, M, offset);
        s2_hat = pskdemod(y2, M, offset);
        
        % Count symbol errors (each symbol error corresponds to b bit errors)
        err_count = err_count + sum(s1_hat ~= s1) + sum(s2_hat ~= s2);
        total_bits = total_bits + length(bits);
        
        %% Record Receiver Constellation for GIF 
        % Concatenate the decoded symbols for visualization
        if packet >= 1 && packet <= 20
            rx_symbols = [y1; y2];
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
    BER(idx) = err_count / total_bits;
end

%% Plot BER vs SNR
figure;
semilogy(SNRdBs, BER, 'b-o', 'LineWidth', 2);
grid on;
xlabel('SNR [dB]');
ylabel('Bit Error Rate (BER)');
axis([SNRdBs([1 end]) 1e-6 1e0]); 
title('BER vs SNR for 2x2 Alamouti Scheme with QPSK');
