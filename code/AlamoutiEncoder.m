function [txSymbols_1, txSymbols_2] = AlamoutiEncoder(X)
% Alamouti STBC Encoder:
% - Inputs:
%   X - modulated symbols matrix (NT x Nsym)
%
% - Outputs:
%   txSymbols_ts1 - Tx symbols for time slot 1 (NT x Nsym)
%   txSymbols_ts2 - Tx symbols for time slot 2 (NT x Nsym)
%
% - Description:
%   Alamouti encoded signal is transmitted from the two transmit
%   antennas over two symbol periods

    NT = size(X,1);
    Nsym = size(X,2);
    if mod(NT,2) ~= 0
        error('Number of transmit antennas must be even.');
    end
    
    txSymbols_1 = zeros(NT, Nsym);
    txSymbols_2 = zeros(NT, Nsym);
    
    % Process each Alamouti pair
    for t = 1:NT/2
        idx1 = 2*t - 1;
        idx2 = 2*t;
        
        % Time slot 1: transmit symbols directly
        txSymbols_1(idx1, :) = X(idx1, :);
        txSymbols_1(idx2, :) = X(idx2, :);
        
        % Time slot 2: apply Alamouti encoding
        txSymbols_2(idx1, :) = -conj(X(idx2, :));
        txSymbols_2(idx2, :) =  conj(X(idx1, :));
    end
end
