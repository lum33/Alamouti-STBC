%% Alamouti_2x1_QPSK_Constellation.m
clear; close all; clc;

%% Parameters
L_frame = 130;          % Number of Alamouti blocks/packet (each block has 2 symbols)
N_Packets = 4000;       % Number of packets
NT = 2; NR = 1;         
M = 4;                  % QPSK 
b = log2(M);            
SNRdBs = 0:2:30;        
offset = pi/4;
txConstDiagram = comm.ConstellationDiagram('Title', 'Transmitter Constellation Diagram');
rxConstDiagram = comm.ConstellationDiagram('Title', 'Receiver Constellation Diagram');
BER = zeros(length(SNRdBs), 1);

%% Simulation loop over SNR
for idx = 1:length(SNRdBs)
    SNRdB = SNRdBs(idx);
    sigma = sqrt(1/(2*10^(SNRdB/10))); % Noise standard deviation
    
    err_count = 0;
    total_bits = 0;
    
    for i_packet = 1:N_Packets
        % Generate random bits for two streams (each stream forms one row of the Alamouti block)
        bits = randi([0 1], 2*L_frame*b, 1);
        % Group bits into symbols (each symbol represented by log2(M) bits)
        data_symbols = bi2de(reshape(bits, b, []).', 'left-msb');
        
        % Split the symbols into two streams
        s1 = data_symbols(1:L_frame);
        s2 = data_symbols(L_frame+1:end);
        
        % QPSK modulation
        x1 = pskmod(s1, M, offset);
        x2 = pskmod(s2, M, offset);
        
        % Independent Rayleigh fading channel for each transmit antenna
        h1 = (randn(L_frame,1) + 1j*randn(L_frame,1))/sqrt(2);
        h2 = (randn(L_frame,1) + 1j*randn(L_frame,1))/sqrt(2);

        % Alamouti encoding: Two time slots per block
        % Time slot 1: antenna 1 sends x1, antenna 2 sends x2.
        r1 = h1.*x1 + h2.*x2;
        % Time slot 2: antenna 1 sends -conj(x2), antenna 2 sends conj(x1).
        r2 = h1.*(-conj(x2)) + h2.*(conj(x1));
        
        % Add complex AWGN noise
        r1 = r1 + sigma*(randn(L_frame,1) + 1j*randn(L_frame,1));
        r2 = r2 + sigma*(randn(L_frame,1) + 1j*randn(L_frame,1));
        
        % Alamouti decoding
        y1 = r1.*conj(h1) + conj(r2).*h2;
        y2 = r1.*conj(h2) - conj(r2).*h1;
        
        % Equalization
        h_power = abs(h1).^2 + abs(h2).^2; % total channel power per block
        y1 = y1 ./ h_power;
        y2 = y2 ./ h_power;
        
        % QPSK demodulation using pskdemod()
        s1_hat = pskdemod(y1, M, offset);
        s2_hat = pskdemod(y2, M, offset);
        
        % Count symbol errors (each symbol error represents an error in log2(M) bits)
        err_count = err_count + sum(s1_hat ~= s1) + sum(s2_hat ~= s2);
        total_bits = total_bits + length(bits);


        tx_symbols = [x1; x2];    % Transmitter constellation points
        rx_symbols = [y1; y2];    % Receiver equalized constellation points
        % Display using the constellation diagram objects
        txConstDiagram(tx_symbols);
        rxConstDiagram(rx_symbols);
    end
    
    % Compute BER for the current SNR
    BER(idx) = err_count / total_bits;
end

%% Plot BER vs SNR
figure;
semilogy(SNRdBs,BER), 
axis([SNRdBs([1 end]) 1e-6 1e0]); 
grid on;  
xlabel('SNR[dB]'); 
ylabel('BER');
title('BER vs SNR for 2x1 Alamouti Scheme with QPSK');
