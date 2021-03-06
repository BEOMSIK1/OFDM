function [H_rv,Y] = SU_MIMO(Data,FFT_Size,Multi_path,GI_Size,Nt,Nr)

Data = randi([0 1],[Nt Data_Size]);
mod_data=base_mod(Data,Modulation_Order);                           % 변조 방식에 따른 데이터 변조
IFFT_data=ifft(mod_data,FFT_Size,2)*sqrt(FFT_Size);                 % 변조된 데이터 IFFT연산(각 반송파에서 보내는 신호의 power를 1로 하기위해 sqrt(반송파 개수)를 곱해줌
Add_CP_data=[IFFT_data(:,FFT_Size-GI_Size+1:end), IFFT_data];       % IFFT연산된 데이터에 CP삽입
copied_cpdata=repmat(Add_CP_data,Nr,1);
h=rayleigh_channel([Multi_path, Nt*Nr]);                             % time domain channel
H=fft(h,FFT_Size,2);                                                % frequency domain channel
for k=1:(Nt*Nr)
    hx(k,:)=conv(copied_cpdata(k,:),h(k,:));
end
for k=1:Nr
    hx_comb(k,:)=sum(hx([(Nt*(k-1)+1):Nt*k],:));
end
y=awgn_noise(hx_comb,SNR(SNR_index));                                    % 채널 통과된 데이터에 awgn추가 (y=h*x+n)
y_remove_CP=y(:,GI_Size+1:GI_Size+FFT_Size);                        %remove CP
Y=fft(y_remove_CP,FFT_Size,2)/sqrt(FFT_Size);                       % CP제거된 데이터 FFT연산(신호의 power를 1로 하기 위해 sqrt(반송파 개수)로 나눠줌)
H_rv=reshape(H,Nt,Nr,[]);                                           % MIMO Channel (transpose 필요)