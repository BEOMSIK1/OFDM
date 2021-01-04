clear,clc, close all

%%모드설정
Mode=[1 2 3 4 5 6]; %line coding 종류
%1:unipolar NRZ,2:Polar NRZ,3:unipolar RZ
%4:polar RZ, 5;bipolar RZ, 6:manchester

Legend_ALL={'Unipolar NRZ','Polar NRZ','Unipolar RZ',...
    'Polar RZ','Bipolar RZ','Manchester'};
Legend_Select=Legend_ALL(Mode);

%%parameter
Rb=1000; %비트율
Tb=1/Rb; %비트 시간
Fs=20000;%샘플링 주파수
Ts=1/Fs; %샘플링 주기

N=Tb/Ts; %한 비트 당 샘플수
%비트 스트림 Data=[11001001];
%L=length[Data]; 비트수
L=1000;%비트수
Data=randi([0 1],[1 L]);%비트스트림
N_ALL=L*N; %number of samples in total bits sector
Index=1:N:N_ALL;%impulse index

t=0:Ts:L*Tb-Ts; %total time
FFT_Size=2^(ceil(log2(length(t))))*2^3; %frequency resolution
d=Fs/FFT_Size;
f=-Fs/2:d:Fs/2-d; %total frequency

%%pulse parameter
ONE=ones(1,N); %NRZ
ONE_ZERO=[ones(1,ceil(N/2)) zeros(1,floor(N/2))]; %RZ
ONE_ONE=[ones(1,ceil(N/2)) -ones(1,floor(N/2))]; %manchester

%%graph setting
t_Ticks=Tb*[0 L];
t_Ticks_Label={'0' 'T_b'};

for k=2 :L
    t_Ticks_Label=[t_Ticks_Label [num2str(k) 'T_b']];
end
%%main
for M=Mode
    s=zeros(1,N_ALL); %시간축 초기화
    Pulse=Data*(rem(M,2)==1)+(2*Data-1)*(rem(M,2)==0);
    if M==5%bipolr rz
        U_=cumsum([1 Pulse]);
        U=(-1).^U_(2:end);
        s(Index)=U.*Pulse;
    else
        s(Index)=Pulse;
    end
    PSF=ONE*(M>=1 && M<=2)+ONE_ZERO*(M>=3&&M<=5)+ONE_ONE*(M==6);
    s=conv(s,PSF);
    s=s(1:N_ALL);
    
    S=fftshift(fft(s,FFT_Size));
    S_Amplitude=abs(S);
    
    figure(M)
    subplot(2,1,1)
    stairs(t,s,'LineWidth',2)
    legend(Legend_Select(M)),set(gca,'FontSize',15),xticks(t_Ticks),xticklabels(t_Ticks_Label)
    subplot(2,1,2)
    plot(f,S_Amplitude,'LineWidth',2)
    legend(Legend_Select(M)),set(gca,'FontSize',15)
    M=M+1;
end




