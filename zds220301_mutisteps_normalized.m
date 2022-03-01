clear all
clc

%%%%%%https://www.zhihu.com/column/control-system ����place����
%FIR fs 16HZ lowpass Fpass = 1Hz  Fstop = 2Hz in 1hz_pass_stop.fcf
%FIR fs 16hz lowpass in 40dbdecay.fcf
h=[0.0005213980539685646                
0.015028694969554298                 
0.030327142272413791                 
0.054429764150456618                 
0.083610269007622362                 
0.11415786162923412                  
0.14118162781628973                  
0.15980849456959786                  
0.16645414020842933                  
0.15980849456959786                  
0.14118162781628973                  
0.11415786162923412                  
0.083610269007622362                 
0.054429764150456618                 
0.030327142272413791                 
0.015028694969554298                 
0.0005213980539685646  ];

% load('fda_30order.mat','NumLen');
% h = NumLen;
len=length(h);
NNUM=50;
frequency=[0.7 ,0.8,0.9,0.95,2];
for index=1:5
for iiii=1:NNUM
iteration_time = 6000;
i=0;
N=len-1;
%initialize
v = zeros(N,iteration_time); w = zeros(N,iteration_time);
x = zeros(N,iteration_time);  %zero state 
e = zeros(1,iteration_time);
filtered_X = zeros(N,iteration_time); 


%Designed Sample Generator
Fs = 16;   % SamplingfrequencyHZ required to be equal to the fs in fdatool                
T = 1/Fs;        % Sampling period  S     
NumSig=iteration_time;     % Number of signal 
t = (0:NumSig-1)*T;        % Time vector
y =sin(2*pi*frequency(index)*t);
%y=sin(2*pi*0.9*t)+ sin(2*pi*0.96*t);
%y=sin(2*pi*1*t);
y = y';

powy=sum(y.^2)/length(y)/Fs

%Noise Generator filtermethod
fs=16;           %��������Ƶ��hz
T_N=iteration_time/fs;%��ʱ��s
t=1/fs:1/fs:T_N;  %ʱ������
L=iteration_time;  %��������=iteration number
sigmaw=powy;          %��������,��λΪdbw
z=wgn(L,1,sigmaw);
[b,a]=butter(8,[0.1/(fs/2),2/(fs/2) ]);%���8�װ�����˹�˲���ϵ����100-200Hz
freqs(b,a)            %���˲�����������
w=filter(b,a,z);

poww=sum(w.^2)/length(w)/Fs

sigmae=2*sigmaw;

[A,B,C,D]=tf2ss(h',ones(1,len));
K=acker(A,B,zeros(1,len-1));

%�ο�����
% sys = ss(A-B*K,B,C,D,1/Fs);
% [y,t,x] = lsim(sys,u,t);

%state and output iteration
for iii=2:iteration_time
   % w = wgn(1,iteration_time,N0*noiseB,'linear');   
    y_cup(iii-1)=C*x(:,iii-1);%+D*e(iii-1);
    x(:,iii)=(A-B*K)*x(:,iii-1)+B*e(iii-1);
    e(iii)=1*D*( y(iii) +w(iii) - C*x(:,iii) )/(sigmaw^2+sigmae^2);   
%%%%%״̬������ôд�����д�����
%     x(:,iii)=A*x(:,iii-1)+B*e(iii-1);
%     e(iii)=D*(y(iii)+w(iii)-C*x(:,iii-1))/(sigmaw^2+sigmae^2);
%     y_cup(iii)=C*x(:,iii-1)+D*e(iii-1); 

end


y_cup=y_cup*8;

 en_w = SampleEntropy(2,0.2*std(w),w);
 en_y = SampleEntropy(2,0.2*std(y_cup),y);
 en_y_cup= SampleEntropy(2,0.2*std(y_cup),y_cup);
 Iyy(index,iiii) = -(en_y_cup - en_w);

end
end

figure(4)
plot(Iyy);title('mutual information');

% for iii=1:iteration_time-1
%     SquDis(iii) = (y_cup(iii)-y(iii))^2;
%     SquDise(iii)= (y_cup(iii)-e(iii))^2;
%     SquDisye(iii)= (e(iii)-y(iii))^2;
% end
 
figure(1)
Y1= fft(y);
%����˫��Ƶ��P2
P1_2 = abs(Y1/NumSig);
%ȡ��ǰ��һ����з���
P1_1 = P1_2(1:NumSig/2+1);
%����ת��Ϊ�����Ƶ
P1_1(2:end-1) = 2*P1_1(2:end-1);
subplot(4,2,1);f = Fs*(0:(NumSig/2))/NumSig;
plot(f,P1_1);title('y��FFT����');
subplot(4,2,2);
plot(y);title('y��ʱ������');


Y_cup=fft(y_cup);
ycupP1_2 = abs(Y_cup/NumSig);
%ȡ��ǰ��һ����з���
ycupP1_1 = ycupP1_2(1:NumSig/2+1);
%����ת��Ϊ�����Ƶ
ycupP1_1(2:end-1) = 2*ycupP1_1(2:end-1);
subplot(4,2,3);f = Fs*(0:(NumSig/2))/NumSig;
plot(f,ycupP1_1);title('ycup��FFT����');
subplot(4,2,4);
plot(y_cup);title('ʵ�����ycup��ʱ������');


E=fft(e);
eP1_2 = abs(E/NumSig);
%ȡ��ǰ��һ����з���
eP1_1 = eP1_2(1:NumSig/2+1);
%����ת��Ϊ�����Ƶ
eP1_1(2:end-1) = 2*eP1_1(2:end-1);
subplot(4,2,5);f = Fs*(0:(NumSig/2))/NumSig;
plot(f,eP1_1);title('e��FFT����');
subplot(4,2,6);
plot(e);title('e��ʱ������');


fftw=fft(w);%����Ҷ�任
P = abs(fftw/NumSig);%ȡ��Ƶ���ԣ�����L
P = P(1:NumSig/2+1);%��ȡǰ���
P(2:end-1)=2*P(2:end-1);%����Ƶ�׷�ֱ�������ǵó���2
subplot(4,2,8)
plot(w)
subplot(4,2,7)
plot(f,P)
title('խ����˹������Ƶ��')

% figure(2)
% subplot(2,2,1);plot(SquDis);title('y��ycup ʱ�� ƽ����');
% subplot(2,2,2);plot(SquDisye);title('y��e ʱ�� ƽ����');
% subplot(2,2,3);plot(SquDise);title('ycup��e ʱ�� ƽ����');

% figure(3)
% waterfall(x);
% 
% figure(5);
% formulation=((C*x+w)*(sigmaw^2/sigmae^2+D^2)-D^2*C*x);
% plot(formulation);




% figure(1)
% waterfall(abs(e)); title('e')
% 
% figure(2)
% waterfall(abs(x)); title('x')
% 
% figure(3)
% waterfall(y); title('y (fixed)')
% 
% figure(4)
% waterfall(abs(filtered_X)); title('filtered x')
