function s_hat = AlamoutiDecoder(Y1, Y2, H)
% Alamouti STBC Decoder:
% Inputs:
%   Y_ts1 - Received symbols at time slot 1 (NR x Nsym)
%   Y_ts2 - Received symbols at time slot 2 (NR x Nsym)
%   H     - Channel frequency response (NR x NT x Nsym)
%
% Outputs:
%   s_hat - Decoded symbols for each transmit antenna (NT x Nsym).
%
% Description:
%   For each sub-block (Alamouti pair) and each symbol time (column) the
%   decoder combines the two received time slots as:
%
%     s1_hat = (sum_rx [ r1*conj(H1) + conj(r2)*H2 ]) / (sum_rx [|H1|^2+|H2|^2])
%     s2_hat = (sum_rx [ r1*conj(H2) - conj(r2)*H1 ]) / (sum_rx [|H1|^2+|H2|^2])

    [NR, Nsym] = size(Y1);
    NT = size(H,2);
    if mod(NT,2) ~= 0
        error('Number of transmit antennas must be even.');
    end
    
    s_hat = zeros(NT, Nsym);
    
    for n = 1:Nsym
        % Process each Alamouti pair
        for t = 1:NT/2
            idx1 = 2*t - 1;
            idx2 = 2*t;
            num1 = 0;
            num2 = 0;
            denom = 0;
            % Combine contributions from all receive antennas
            for rx = 1:NR
                H1 = H(rx, idx1, n);
                H2 = H(rx, idx2, n);
                r1 = Y1(rx, n);
                r2 = Y2(rx, n);
                num1 = num1 + r1*conj(H1) + conj(r2)*H2;
                num2 = num2 + r1*conj(H2) - conj(r2)*H1;
                denom = denom + abs(H1)^2 + abs(H2)^2;
            end
            % Decode the pair of symbols
            s_hat(idx1, n) = num1 / denom;
            s_hat(idx2, n) = num2 / denom;
        end
    end
end
