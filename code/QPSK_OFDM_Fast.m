%% Parameter Setup
L_frame = 1024;      % Increase total QPSK symbols per packet to get more OFDM symbols
N_Packets = 10000;   % Number of packets remains the same
M = 4;               % QPSK modulation order
b = log2(M);         % Bits per symbol
offset = pi/4;       % Phase offset for QPSK symbols  
SNRdBs = -5:2:25;    % Range of SNR values (dB)
total_bits_qpsk = L_frame * b; % Total bits per packet

% OFDM Parameters 
N_subcarriers = 64;  % Keep subcarriers the same
CP_length = 16;      % Reduce CP length since channel is flat fading
N_OFDM = L_frame / N_subcarriers;  % Now equals 16 OFDM symbols per packet

figure(1)
%% Rayleigh Channel Configuration
rayleighChan = comm.RayleighChannel(...
    'SampleRate', 1e6, ...                 
    'PathDelays', [0 1e-6 3e-6 5e-6], ...   
    'AveragePathGains', [0 -3 -6 -9], ...   
    'MaximumDopplerShift', 100, ...          
    'PathGainsOutputPort', true);            

%{
rayleighChan = comm.RayleighChannel(...
    'SampleRate', 1e6, ...             % Set sample rate (adjust as needed)
    'PathDelays', 0, ...          % Two paths: one at 0 and one with a delay of 1 microsecond
    'AveragePathGains', 0, ...      % Average gains (in dB) for the two paths
    'MaximumDopplerShift', 100, ...      % Fast fading: Doppler shift in Hz (adjust for different mobility scenarios)
    'PathGainsOutputPort', true);        % Enable output of individual path gains

%}
%% QPSK Simulation with Fast Fading Multipath Channel
figure(1)
BER_qpsk = zeros(length(SNRdBs), 1);

for idx = 1:length(SNRdBs)
    SNRdB = SNRdBs(idx);
    sigma = sqrt(1/(2*10^(SNRdB/10))); % AWGN noise standard deviation
    err_count = 0;
    
    % Reset the channel object at the start of each SNR point to clear its state
    reset(rayleighChan);
    
    for packet = 1:N_Packets
        % Generate random data bits and map to QPSK symbols
        bits = randi([0 1], total_bits_qpsk, 1);
        X_int = bi2de(reshape(bits, b, []).', 'left-msb');
        X_mod = pskmod(X_int, M, offset);
        
        % Pass the modulated symbols through the Rayleigh multipath channel.
        % The channel object outputs both the channel-affected signal and the individual path gains.
        [Y_channel, pathGains] = rayleighChan(X_mod);  % Note: X_mod is a column vector here
        Y_channel = Y_channel.';  % Convert to row vector for further processing
        
        % Add AWGN noise
        noise = sigma * (randn(1, L_frame) + 1j*randn(1, L_frame));
        Y_noisy = Y_channel + noise;
        
        % For coherent detection, assume perfect channel knowledge.
        % Here, we use the gain of the first tap for equalization.
        % (In a more advanced receiver, you might combine taps or use channel estimation.)
        % Since 'pathGains' is an L_frame-by-2 matrix (one column per path) and was output
        % corresponding to the input symbols, we take the first column.
        pathGains = pathGains(:,1).'; % Convert first tap gains to a row vector
        Y_eq = Y_noisy ./ pathGains;
        
        % Maximum Likelihood Detection for QPSK symbols
        refSymbols = pskmod((0:M-1).', M, offset);
        s_hat_ml = zeros(L_frame,1);
        for n = 1:L_frame
            y = Y_eq(n);
            distances = abs(y - refSymbols).^2;
            [~, minIdx] = min(distances);
            s_hat_ml(n) = refSymbols(minIdx);
        end
        
        % Demodulation to bits and error counting
        X_hat = pskdemod(s_hat_ml, M, offset);
        err_count = err_count + sum(X_hat ~= X_int);
    end
    
    BER_qpsk(idx) = err_count / (total_bits_qpsk * N_Packets);
end

% Plot the Bit Error Rate (BER) performance

semilogy(SNRdBs, BER_qpsk, 'g-^', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
title('BER Performance of QPSK with severe Multipath Fading Channel');
hold on;

%% 2. QPSK Simulation with OFDM (Block Fading per OFDM Symbol)
BER_qpsk_OFDM = zeros(length(SNRdBs), 1);

for idx = 1:length(SNRdBs)
    SNRdB = SNRdBs(idx);
    sigma = sqrt(1/(2*10^(SNRdB/10))); % Noise standard deviation
    err_count = 0;
    
    % Reset the channel object for each SNR point
    reset(rayleighChan);
    
    for packet = 1:N_Packets
        % Generate random bits and QPSK modulation for one OFDM packet
        bits = randi([0 1], total_bits_qpsk, 1);
        X_int = bi2de(reshape(bits, b, []).', 'left-msb');
        X_mod = pskmod(X_int, M, offset);
        % Reshape modulated symbols into OFDM blocks
        X_mod_matrix = reshape(X_mod, N_subcarriers, N_OFDM);
        
        % Initialize containers for equalized symbols per packet
        Y_eq_OFDM = zeros(N_subcarriers, N_OFDM);
        
        % Process each OFDM symbol individually (block fading assumption)
        for k = 1:N_OFDM
            % --- OFDM Modulation for symbol k ---
            ifft_data = ifft(X_mod_matrix(:, k), N_subcarriers);
            tx_symbol = [ifft_data(end-CP_length+1:end); ifft_data];
            tx_symbol = tx_symbol / sqrt(mean(abs(tx_symbol).^2)); % Normalize symbol power
        
            % --- Channel: Pass the symbol through the Rayleigh channel ---
            [y_sym, pathGains_sym] = rayleighChan(tx_symbol);
            
            % --- Compute frequency response per OFDM symbol ---
            % Here we use the channel impulse response from the first sample as an example.
            h_impulse = pathGains_sym(1, :).';  % Extract channel gains (one per path)
            if length(h_impulse) < N_subcarriers
                h_impulse = [h_impulse; zeros(N_subcarriers - length(h_impulse), 1)];
            end
            H_sym = fft(h_impulse, N_subcarriers);
            
            % --- Add AWGN noise ---
            noise_sym = sigma * (randn(size(y_sym)) + 1j*randn(size(y_sym)));
            y_sym_noisy = y_sym + noise_sym;
            
            % --- OFDM Demodulation for symbol k ---
            rx_noCP = y_sym_noisy(CP_length+1:end);
            Y_sym = fft(rx_noCP, N_subcarriers);
            
            % --- Equalize per subcarrier using the frequency response ---
            Y_eq_OFDM(:, k) = Y_sym ./ H_sym;
        end

        
        % Serialize the equalized OFDM symbols into a vector for detection
        Y_eq_vec = Y_eq_OFDM(:);
        
        % --- QPSK Detection ---
        refSymbols = pskmod((0:M-1).', M, offset);
        s_hat_ml = zeros(length(Y_eq_vec), 1);
        for n = 1:length(Y_eq_vec)
            y = Y_eq_vec(n);
            distances = abs(y - refSymbols).^2;
            [~, minIdx] = min(distances);
            s_hat_ml(n) = refSymbols(minIdx);
        end
        
        % Demodulation to bits and error counting
        X_hat = pskdemod(s_hat_ml, M, offset);
        err_count = err_count + sum(X_hat ~= X_int);
    end
    BER_qpsk_OFDM(idx) = err_count / (total_bits_qpsk * N_Packets);
end

semilogy(SNRdBs, BER_qpsk_OFDM, 'b-o', 'LineWidth', 2);
legend('non-OFDM','OFDM');
saveas(figure(1), 'QPSK_OFDM_Fast2.jpg')