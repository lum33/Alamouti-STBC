function [tmp1, sym_tab, P] = modulator(bits, b)
    M = 2^b; % Modulation order (e.g., QPSK: M = 4)
    sym_tab = exp(1j * 2 * pi * (0:M-1) / M); % Symbol table (constellation points)
    P = M; % Number of constellation points
    indices = bi2de(reshape(bits, log2(M), []).') + 1; % Convert bits to symbols
    tmp1 = sym_tab(indices); % Map bits to symbols
end
