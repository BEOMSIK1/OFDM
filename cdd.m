clc, clear, close all;
%% Parameters

Modulation_Order=2;                   % ���� ��� 2:QPSK, 4:QAM
FFT_Size=128;                         % �ݼ��� ����
Data_Size =FFT_Size*Modulation_Order; % ������ ũ��
GI_Size=FFT_Size/4;                   % CP size
Multi_path=7;                         % ���߰�� ����
N=2;                                  % �۽� ���׳� ����
M=1;                                  % ���� ���׳� ����
cl=3;                                 % constraint length
code_rate=1/2;                        % code rate
                             
SNR=0:3:30;
Iteration=5000;

%% CDD
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        %% MISO-OFDM (CDD�� ���� �񱳴��)
        Data = randi([0 1],[N Data_Size]);     
        mod_data=base_mod(Data,Modulation_Order);                           % ���� ��Ŀ� ���� ������ ����
        IFFT_data=ifft(mod_data,FFT_Size,2)*sqrt(FFT_Size);                 % ������ ������ IFFT����(�� �ݼ��Ŀ��� ������ ��ȣ�� power�� 1�� �ϱ����� sqrt(�ݼ��� ����)�� ������
        Add_CP_data=[IFFT_data(:,FFT_Size-GI_Size+1:end), IFFT_data];       % IFFT����� �����Ϳ� CP����
 
        h=rayleigh_channel([Multi_path N]);                                 % time domain channel
        H=fft(h,FFT_Size,2);                                                % frequency domain channel
        
        for Num_channel=1:N
            hx(Num_channel,:)=conv(Add_CP_data(Num_channel,:),h(Num_channel,:));  % ���۵� �����Ͱ� ���߰�� ä�� ���(h*x)
        end
        
        y=awgn_noise(hx,SNR(SNR_index));                                    % ä�� ����� �����Ϳ� awgn�߰� (y=h*x+n)
        y_remove_CP=y(:,GI_Size+1:GI_Size+FFT_Size);                        %remove CP
        Y=fft(y_remove_CP,FFT_Size,2)/sqrt(FFT_Size);                       % CP���ŵ� ������ FFT����(��ȣ�� power�� 1�� �ϱ� ���� sqrt(�ݼ��� ����)�� ������)
        Y_equalize=Y./H; 
        Y_demod=base_demod(Y_equalize,Modulation_Order);                 % equalize�� �����͸� ���� 
            


        num_error(Iter,SNR_index)=biterr(Data,Y_demod);                     % �� SNR�� ��Ʈ���������� Iteration���� �迭�� ����
    end
end
error_rate=(sum(num_error,1)/(Data_Size*N))/Iteration;

semilogy(SNR,error_rate,'-o')
title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
grid on