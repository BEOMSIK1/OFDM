clc, clear, close all;

%% Parameters
FFT_Size=128;                             % �ݼ��� ����
GI_Size=FFT_Size/4;                       % CP size
Modulation_Order=2;                       % 1:BPSK  2:QPSK  4: 16QAM  6: 64QAM  8: 256QAM
Data_Size=FFT_Size*Modulation_Order; 
Multi_path=7;                             % ���߰�� ����
Nt=2;                                     % �۽ž��׳� ���� 
Nr=2;                                     % ���ž��׳� ����
SNR=0:3:30;
Iteration=1000;
%% MIMO OFDM System
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        Data=randi([0 1],[Nt Data_Size]);                             % 0 or 1�� ���� ������ Data_Sizeũ�⸸ŭ�� ������ ����
        mod_data=base_mod(Data,Modulation_Order);                     % ���� ��Ŀ� ���� ������ ����
        IFFT_data=ifft(mod_data)*sqrt(FFT_Size);                      % ������ ������ IFFT����(�� �ݼ��Ŀ��� ������ ��ȣ�� power�� 1�� �ϱ����� sqrt(�ݼ��� ����)�� ������
        Add_CP_data=[IFFT_data(:,FFT_Size-GI_Size+1:end), IFFT_data]; % IFFT����� �����Ϳ� CP����

        h=rayleigh_channel(Multi_path);                               % ���߰�� ä��
        H=fft(h,FFT_Size);                                            % ä�� ���ļ� ����
        for K=1:Nt
            hx(K,:)=conv(Add_CP_data(K,:),h);                         % ���۵� �����Ͱ� ���߰�� ä���� ���
        
            y(K,:)=awgn_noise(hx(K,:),SNR(SNR_index));                              % ä�� ����� �����Ϳ� awgn�߰�
        
            y_remove_CP(K,:)=y(K,GI_Size+1:GI_Size+FFT_Size);                  % remove CP
            Y(K,:)=fft(y_remove_CP(K,:),FFT_Size)/sqrt(FFT_Size);                  % CP���ŵ� ������ FFT����(��ȣ�� power�� 1�� �ϱ� ���� sqrt(�ݼ��� ����)�� ������)
        
            Y_equalize(K,:)=Y(K,:)./H;                                              % ��ȭ����
            Y_demod(K,:)=base_demod(Y_equalize(K,:),Modulation_Order);              % equalize�� �����͸� ���� 
        end
        
        
        num_error(Iter,SNR_index)=biterr(Data,Y_demod);               % �� SNR�� ��Ʈ���������� Iteration���� �迭�� ����
    end
 
end
error_rate=(sum(num_error,1)/(Data_Size*Nt))/Iteration;
%% graph
semilogy(SNR,error_rate,'-o')
title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
grid on