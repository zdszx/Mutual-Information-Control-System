clear all
clc

%%%%%%https://www.zhihu.com/column/control-system 采用place函数
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
fs=16;           %噪声采样频率hz
T_N=iteration_time/fs;%总时间s
t=1/fs:1/fs:T_N;  %时间向量
L=iteration_time;  %样本数量=iteration number
sigmaw=powy;          %噪声功率,单位为dbw
z=wgn(L,1,sigmaw);
[b,a]=butter(8,[0.1/(fs/2),2/(fs/2) ]);%获得8阶巴特沃斯滤波器系数，100-200Hz
freqs(b,a)            %画滤波器特性曲线
w=filter(b,a,z);

poww=sum(w.^2)/length(w)/Fs

sigmae=2*sigmaw;

[A,B,C,D]=tf2ss(h',ones(1,len));
K=acker(A,B,zeros(1,len-1));

%参考函数
% sys = ss(A-B*K,B,C,D,1/Fs);
% [y,t,x] = lsim(sys,u,t);

%state and output iteration
for iii=2:iteration_time
   % w = wgn(1,iteration_time,N0*noiseB,'linear');   
    y_cup(iii-1)=C*x(:,iii-1);%+D*e(iii-1);
    x(:,iii)=(A-B*K)*x(:,iii-1)+B*e(iii-1);
    e(iii)=1*D*( y(iii) +w(iii) - C*x(:,iii) )/(sigmaw^2+sigmae^2);   
%%%%%状态方程怎么写不会有大区别
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
%计算双侧频谱P2
P1_2 = abs(Y1/NumSig);
%取出前面一半进行分析
P1_1 = P1_2(1:NumSig/2+1);
%最终转化为单侧幅频
P1_1(2:end-1) = 2*P1_1(2:end-1);
subplot(4,2,1);f = Fs*(0:(NumSig/2))/NumSig;
plot(f,P1_1);title('y的FFT分析');
subplot(4,2,2);
plot(y);title('y的时域序列');


Y_cup=fft(y_cup);
ycupP1_2 = abs(Y_cup/NumSig);
%取出前面一半进行分析
ycupP1_1 = ycupP1_2(1:NumSig/2+1);
%最终转化为单侧幅频
ycupP1_1(2:end-1) = 2*ycupP1_1(2:end-1);
subplot(4,2,3);f = Fs*(0:(NumSig/2))/NumSig;
plot(f,ycupP1_1);title('ycup的FFT分析');
subplot(4,2,4);
plot(y_cup);title('实际输出ycup的时域序列');


E=fft(e);
eP1_2 = abs(E/NumSig);
%取出前面一半进行分析
eP1_1 = eP1_2(1:NumSig/2+1);
%最终转化为单侧幅频
eP1_1(2:end-1) = 2*eP1_1(2:end-1);
subplot(4,2,5);f = Fs*(0:(NumSig/2))/NumSig;
plot(f,eP1_1);title('e的FFT分析');
subplot(4,2,6);
plot(e);title('e的时域序列');


fftw=fft(w);%傅里叶变换
P = abs(fftw/NumSig);%取幅频特性，除以L
P = P(1:NumSig/2+1);%截取前半段
P(2:end-1)=2*P(2:end-1);%单侧频谱非直流分量记得乘以2
subplot(4,2,8)
plot(w)
subplot(4,2,7)
plot(f,P)
title('窄带高斯噪声（频域）')

% figure(2)
% subplot(2,2,1);plot(SquDis);title('y与ycup 时序 平方差');
% subplot(2,2,2);plot(SquDisye);title('y与e 时序 平方差');
% subplot(2,2,3);plot(SquDise);title('ycup与e 时序 平方差');

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
