clc, clear, close all;
%% Parameters
% Nt = 100;                                   % 송신안테나 개수
% Nr = 1;                                     % 수신안테나 개수 (유저 당)
% K = 10;                                     % 유저의 수
modulation_order = 4;                       % 1:BPSK  2:QPSK  4: 16QAM  6: 64QAM  8: 256QAM
symbol_size = 128;
data_size = symbol_size*modulation_order;
P=-10:3:20;                                   % power
SNR=0;
iteration=1000;
%% Massive MIMO
for power_index=1:length(P)
    for iter=1:iteration
        %% Massive MIMO
        [num_error_100x10,S1] = Massive_MIMO_ZF(P,power_index,100,10,SNR,modulation_order,data_size);
        num_error_ZF_100x10(iter,power_index)=num_error_100x10;
        [num_error_100x20,S2] = Massive_MIMO_ZF(P,power_index,100,20,SNR,modulation_order,data_size);
        num_error_ZF_100x20(iter,power_index)=num_error_100x20;
        [num_error_200x20] = Massive_MIMO_ZF(P,power_index,200,20,SNR,modulation_order,data_size);
        num_error_ZF_200x20(iter,power_index)=num_error_200x20;
        [num_error_200x40] = Massive_MIMO_ZF(P,power_index,200,40,SNR,modulation_order,data_size);
        num_error_ZF_200x40(iter,power_index)=num_error_200x40;
        %% Sum rate
        SR_100x10(iter,power_index)=S1;
        SR_100x20(iter,power_index)=S2;
        
        

    end
end
error_rate_ZF_100x10=(sum(num_error_ZF_100x10,1)/(data_size*10))/iteration;
error_rate_ZF_100x20=(sum(num_error_ZF_100x20,1)/(data_size*20))/iteration;
error_rate_ZF_200x20=(sum(num_error_ZF_200x20,1)/(data_size*20))/iteration;
error_rate_ZF_200x40=(sum(num_error_ZF_200x40,1)/(data_size*40))/iteration;

sum_rate_100x10=(sum(abs(SR_100x10),1))/iteration;
sum_rate_100x20=(sum(abs(SR_100x20),1))/iteration;
%% graph
semilogy(P,error_rate_ZF_100x10,'-o')
hold on
semilogy(P,error_rate_ZF_100x20,'-o')
hold on
semilogy(P,error_rate_ZF_200x20,'-o')
hold on
semilogy(P,error_rate_ZF_200x40,'-o')
hold on
title('BER Performance'), xlabel('Transmit Power(dB)'),ylabel('BER'),legend('100x10 (16-QAM)','100x20 (16-QAM)','200x20 (16-QAM)','200x40 (16-QAM)')
grid on

figure
plot(P,sum_rate_100x10,'-o')
hold on
plot(P,sum_rate_100x20,'-o')
title('Sum Rate'), xlabel('Transmit Power(dB)'),ylabel('Sum Rate(bps/Hz)'),legend('100x10','100x20')
grid on
