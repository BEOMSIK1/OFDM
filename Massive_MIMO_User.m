clc, clear, close all;
%% Parameters
% Na                                    % 송신안테나 개수
% Nr=1                                  % 수신안테나 개수 (유저 당)
% K                                     % 유저의 수
modulation_order = 4;                       % 1:BPSK  2:QPSK  4: 16QAM  6: 64QAM  8: 256QAM
symbol_size = 128;
data_size = symbol_size*modulation_order;
P = [5 10];                                % fixed power (5,10dB) / power_index=1 -> 5dB / power_index=2 -> 10dB 
Na = 200;                                  % Tx antenna     
K = [1 5 10 20 : 20 : 200];
SNR=0;
iteration=1000;
%% Massive MIMO
for user_index=1:length(K)
    for iter=1:iteration
        %% Massive MIMO
        [num_error_200xK_5db,S1] = Massive_MIMO_ZF(P,1,Na,K(user_index),SNR,modulation_order,data_size);
        num_error_ZF_200xK_5db(iter,user_index)=num_error_200xK_5db./K(user_index);
        [num_error_200xK_10db,S2] = Massive_MIMO_ZF(P,2,Na,K(user_index),SNR,modulation_order,data_size);
        num_error_ZF_200xK_10db(iter,user_index)=num_error_200xK_10db./K(user_index);
        %% Sum rate
        SR_200xK_5db(iter,user_index)=S1;
        SR_200xK_10db(iter,user_index)=S2;
        
        

    end
end
error_rate_ZF_200xK_5db=(sum(num_error_ZF_200xK_5db,1)/(data_size))/iteration;
error_rate_ZF_200xK_10db=(sum(num_error_ZF_200xK_10db,1)/(data_size))/iteration;

sum_rate_200xK_5db=(sum(abs(SR_200xK_5db),1))/iteration;
sum_rate_200xK_10db=(sum(abs(SR_200xK_10db),1))/iteration;
%% graph
semilogy(K,error_rate_ZF_200xK_5db,'-o')
hold on
semilogy(K,error_rate_ZF_200xK_10db,'-o')
title('BER Performance'), xlabel('Number of user'),ylabel('BER'),legend('Nt=200 P=5dB (16-QAM)','Nt=200 P=10dB (16-QAM)')
grid on

figure
plot(K,sum_rate_200xK_5db,'-o')
hold on
plot(K,sum_rate_200xK_10db,'-o')
title('Sum Rate'), xlabel('Number of user'),ylabel('Sum Rate(bps/Hz)'),legend('Nt=200 P=5dB','Nt=200 P=10dB')
grid on
