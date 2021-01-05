clc, clear, close all

%% parameter
f0=1000; %Fundamental Frequency
T=1/f0; %Symbol Duration
L=10; %Number of Data
Data=randi([0 1], [1 L]); %Random Bits
Index=find(Data==0); %Index of Bit 0
Data(Index)=-1; %Change 0 to -1
t=linspace(0,T,1000);

%%OFDM signal implementation
OFDM_Cos=0;
OFDM_Sin=0; %OFDM signal initialization
for k=1:L
   Cos_carrier=Data(k)*cos(2*pi*k*f0*t);
   Sin_carrier=Data(k)*sin(2*pi*k*f0*t);
   OFDM_Cos=OFDM_Cos+Cos_carrier;
   OFDM_Sin=OFDM_Sin+Sin_carrier;
end

N_IDFT=128;
N_IFFT=128;

%%IDFT implementation
IDFT_sum=0;
for n=1:N_IDFT
    for k=1:N_IDFT
        Data=[Data zeros(1,N_IDFT-length(Data))];
        IDFT_sum(k)=Data(k)*exp((j*2*pi*k*n)/N_IDFT);
        IDFT(n)=sum(IDFT_sum);
    end
end

%%plot IDFT and OFDM graph
figure
stem(linspace(0,T,N_IDFT),IDFT)
hold on
plot(t,OFDM_Cos)
xlabel('t [sec]'), ylabel('y(t)'), grid on,legend('IDFT Cos','OFDM Cos')


figure
stem(linspace(0,T,N_IDFT),imag(IDFT))
hold on
plot(t,OFDM_Sin)
xlabel('t [sec]'), ylabel('y(t)'), grid on,legend('IDFT Sin','OFDM Sin')

%%IFFT implementation
IFFT_sum=0;
for n=1:N_IFFT
   for k=1:N_IFFT
       Data=[Data(1:L) zeros(1,N_IFFT-length(Data(1:L)))];
       IFFT_sum(k)=Data(k)*exp((j*2*pi*k*n)/N_IFFT);
       IFFT(n)=sum(IFFT_sum);
   end
end

%IFFT=ifft(Data(1:L),N_IFFT)*N_IFFT;

%%plot IFFT and ODFM graph
figure
stem(linspace(0,T,N_IFFT),IFFT)
hold on
plot(t,OFDM_Cos)
xlabel('t [sec]'), ylabel('y(t)'), grid on,legend('IFFT Cos','OFDM Cos')

figure
stem(linspace(0,T,N_IFFT),imag(IFFT))
hold on
plot(t,OFDM_Sin)
xlabel('t [sec]'), ylabel('y(t)'), grid on,legend('IFFT Sin','OFDM Sin')