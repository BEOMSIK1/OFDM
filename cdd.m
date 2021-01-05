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
Iteration=1000;

%% CDD
for SNR_index=1:length(SNR)
    for Iter=1:Iteration
        %% MISO-OFDM (CDD�� ���� �񱳴��)
        Data = randi([0 1],[N Data_Size]);     
        mod_data=base_mod(Data,Modulation_Order);                           % ���� ��Ŀ� ���� ������ ����
        IFFT_data=ifft(mod_data)*sqrt(FFT_Size);                            % ������ ������ IFFT����(�� �ݼ��Ŀ��� ������ ��ȣ�� power�� 1�� �ϱ����� sqrt(�ݼ��� ����)�� ������
        Add_CP_data=[IFFT_data(:,FFT_Size-GI_Size+1:end), IFFT_data];       % IFFT����� �����Ϳ� CP����

        for Num_channel=1:N                                                 % �۽� ���׳� ������ŭ ä�� ����
            
            h(Num_channel,:)=rayleigh_channel(Multi_path);          % time domain channel
            H(Num_channel,:)=fft(h(Num_channel,:),FFT_Size);        % frequency domain channel
            
            hx(Num_channel,:)=conv(Add_CP_data(Num_channel,:),h(Num_channel,:));  % ���۵� �����Ͱ� ���߰�� ä�� ���(h*x)
            y(Num_channel,:)=awgn_noise(hx(Num_channel,:),SNR(SNR_index)); % ä�� ����� �����Ϳ� awgn�߰� (y=h*x+n)
            
            y_remove_CP(Num_channel,:)=y(Num_channel,GI_Size+1:GI_Size+FFT_Size);      %remove CP
            Y(Num_channel,:)=fft(y_remove_CP(Num_channel,:),FFT_Size)/sqrt(FFT_Size);  % CP���ŵ� ������ FFT����(��ȣ�� power�� 1�� �ϱ� ���� sqrt(�ݼ��� ����)�� ������)
            Y_equalize(Num_channel,:)=Y(Num_channel,:)./H(Num_channel,:); 
            Y_demod(Num_channel,:)=base_demod(Y_equalize(Num_channel,:),Modulation_Order);                    % equalize�� �����͸� ���� 
            

        end
            num_error(Iter,SNR_index)=biterr(Data,Y_demod);                     % �� SNR�� ��Ʈ���������� Iteration���� �迭�� ����
    end
end