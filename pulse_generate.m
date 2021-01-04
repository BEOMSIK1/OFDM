clear,clc,close all
%%parameter in time
Fs=10; %���ø� ���ļ�
Ts=1/Fs; %���ø��ֱ�

W=5;%�ð�����
t=-W:Ts:W-Ts; %��ü �ð�

%%����A, �� 2T�� �簢�޽�
A=1;
T=0.5;
s=zeros(1,length(t));
Index=find(t>=-T&t<T);
s(Index)=A;

%% parameter in frequency

FFT_Size=2048; %���ļ� �ػ�
d=Fs/FFT_Size;
f=-Fs/2:d:Fs/2-d; %��ü ���ļ�

%%����Ʈ��

S=fftshift(fft(s,FFT_Size)*Ts);
S_Amplitude=abs(S);

%%�׷���
t_Ticks=[-T T];
t_Ticks_Label={'-T' 'T'};
f_Ticks=-Fs/2:1/(2*T):Fs/2;
f_Minus_Ticks_Label={};
f_Plus_Ticks_Label={};

for k=1:length(1/(2*T):1/(2*T):Fs/2)
    f_Minus_Ticks_label=[[num2str(-k)'/T'] f_Minus_Ticks_Label];
    f_Plus_Ticks_Label=[f_Plus_Ticks_Label [num2str(k) '/T']];
end
f_Ticks_Label=[f_Minus_Ticks_Label '0' f_Plus_Ticks_Label];

figure(1)
stairs(t,s,'LineWidth',2)
xticks(t_Ticks),xticklabels(t_Ticks_Label),grid on, set(gca,'FontSize',15)
figure(3)
plot(f,S_Amplitude,'LineWidth',2)
grid on, set(gca,'FontSize',15),xticks(f_Ticks),xticklabels(f_Ticks_Label)
