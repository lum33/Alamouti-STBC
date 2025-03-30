% Alamouti_scheme.m
clear; clf;
%% Initialization
L_frame=130; % Number of frames/packet
N_Packets=4000;  % Number of packets 
NT=2; 
NR=1; 
b=2; % Number of bits per symbol = QPSK
SNRdBs=0:2:30; 
sq_NT=sqrt(NT); 
sq2=sqrt(2);
% Simulation loop over SNR 
for i_SNR=1:length(SNRdBs)
    SNRdB=SNRdBs(i_SNR);  
    sigma=sqrt(0.5/(10^(SNRdB/10))); % Noise standard deviation

    % Packet transmission simulation
    for i_packet=1:N_Packets
        % randi([minVal, maxVal], m, n) -> m x n matrix
        msg_symbol=randi([0, 1], L_frame*b,NT); 
        tx_bits=msg_symbol.';  
        tmp=[];   
        tmp1=[];

        % Modulation
        for i=1:NT
          [tmp1,sym_tab,P] = modulator(tx_bits(i,:),b); ;
          tmp=[tmp; tmp1];
        end
        
        % Alamouti encoding
        X=tmp.'; 
        X1=X; 
        X2=[-conj(X(:,2)) conj(X(:,1))];

        % Channel simulation
        for n=1:NT
            Hr(n,:,:)=(randn(L_frame,NT)+j*randn(L_frame,NT))/sq2;
        end
        H=reshape(Hr(n,:,:),L_frame,NT);  
        Habs(:,n)=sum(abs(H).^2,2);

        % Received signal
        R1 = sum(H.*X1,2)/sq_NT + sigma*(randn(L_frame,1)+j*randn(L_frame,1));
        R2 = sum(H.*X2,2)/sq_NT + sigma*(randn(L_frame,1)+j*randn(L_frame,1));

        % Alamouti decoding
        Z1 = R1.*conj(H(:,1)) + conj(R2).*H(:,2);
        Z2 = R1.*conj(H(:,2)) - conj(R2).*H(:,1);

        % Demodulation and decision
        for m=1:P
            tmp = (-1+sum(Habs,2))*abs(sym_tab(m))^2;
            d1(:,m) = abs(sum(Z1,2)-sym_tab(m)).^2 + tmp;
            d2(:,m) = abs(sum(Z2,2)-sym_tab(m)).^2 + tmp;
        end

        % BER
        [y1,i1]=min(d1,[],2);   S1d=sym_tab(i1).';    clear d1
        [y2,i2]=min(d2,[],2);   S2d=sym_tab(i2).';    clear d2
        Xd = [S1d S2d];   tmp1=X>0;  tmp2=Xd>0;
        noeb_p(i_packet) = sum(sum(tmp1~=tmp2));% for coded
    
    end % End of FOR loop for i_packet
    BER(i_SNR) = sum(noeb_p)/(N_Packets*L_frame*b);
end    % End of FOR loop for i_SNR

% Plot
semilogy(SNRdBs,BER), 
axis([SNRdBs([1 end]) 1e-6 1e0]); 
grid on;  
xlabel('SNR[dB]'); 
ylabel('BER');