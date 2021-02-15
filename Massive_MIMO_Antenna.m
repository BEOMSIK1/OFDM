clc, clear, close all;
%% Parameters
% Na                                    % 송신안테나 개수
% Nr=1                                  % 수신안테나 개수 (유저 당)
% K                                     % 유저의 수
modulation_order = 4;                       % 1:BPSK  2:QPSK  4: 16QAM  6: 64QAM  8: 256QAM
symbol_size = 128;
data_size = symbol_size*modulation_order;
P = 10;                                % fixed power (10dB)
Na = 20:20:200;                        % Tx antenna                 
SNR=0;
iteration=1000;
%% Massive MIMO
for antenna_index=1:length(Na)
    power_index = 1;
    for iter=1:iteration
        %% Massive MIMO
        [num_error_NaX10,S1] = Massive_MIMO_ZF(P,power_index,Na(antenna_index),10,SNR,modulation_order,data_size);
        num_error_ZF_NaX10(iter,antenna_index)=num_error_NaX10;
        [num_error_NaX20,S2] = Massive_MIMO_ZF(P,power_index,Na(antenna_index),20,SNR,modulation_order,data_size);
        num_error_ZF_NaX20(iter,antenna_index)=num_error_NaX20;
        %% Sum rate
        SR_NaX10(iter,antenna_index)=S1;
        SR_NaX20(iter,antenna_index)=S2;
        
        

    end
end
error_rate_ZF_NaX10=(sum(num_error_ZF_NaX10,1)/(data_size*10))/iteration;
error_rate_ZF_NaX20=(sum(num_error_ZF_NaX20,1)/(data_size*20))/iteration;

sum_rate_NaX10=(sum(abs(SR_NaX10),1))/iteration;
sum_rate_NaX20=(sum(abs(SR_NaX20),1))/iteration;
%% graph
semilogy(Na,error_rate_ZF_NaX10,'-o')
hold on
semilogy(Na,error_rate_ZF_NaX20,'-o')
title('BER Performance'), xlabel('Number of Tx antenna'),ylabel('BER'),legend('K=10 P=10dB (16-QAM)','K=20 P=10dB (16-QAM)')
grid on

figure
plot(Na,sum_rate_NaX10,'-o')
hold on
plot(Na,sum_rate_NaX20,'-o')
title('Sum Rate'), xlabel('Number of Tx antenna'),ylabel('Sum Rate(bps/Hz)'),legend('K=10','K=20')
grid on
