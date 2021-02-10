function [Y] = MISO_OFDM(mod_data,FFT_Size,GI_Size,h,SNR_index,Nt)

SNR=0:3:30;
IFFT_data=ifft(mod_data,FFT_Size,2)*sqrt(FFT_Size)/sqrt(Nt);                 % 변조된 데이터 IFFT연산(각 반송파에서 보내는 신호의 power를 1로 하기위해 sqrt(반송파 개수)를 곱해줌
Add_CP_data=[IFFT_data(:,FFT_Size-GI_Size+1:end), IFFT_data];       % IFFT연산된 데이터에 CP삽입


for Num_channel=1:Nt
    hx(Num_channel,:)=conv(Add_CP_data(Num_channel,:),h(Num_channel,:));  % 전송된 데이터가 다중경로 채널 통과(h*x)
end

y=awgn_noise(hx,SNR(SNR_index));                                    % 채널 통과된 데이터에 awgn추가 (y=h*x+n)
y_remove_CP=y(:,GI_Size+1:GI_Size+FFT_Size);                        %remove CP
Y=fft(y_remove_CP,FFT_Size,2)/sqrt(FFT_Size)*sqrt(Nt);                       % CP제거된 데이터 FFT연산(신호의 power를 1로 하기 위해 sqrt(반송파 개수)로 나눠줌)