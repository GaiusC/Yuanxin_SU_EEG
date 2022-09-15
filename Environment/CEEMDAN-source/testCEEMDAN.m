%CEEMDAN方法的测试脚本程序
%% 1.生成仿真信号
fs = 400;  %采样频率
t = 0:1/fs:0.75; %时间轴
x = sin(2*pi*4*t); %低频正弦信号
y = 0.5*sin(2*pi*120*t); %高频正弦信号
for i = 1:length(t) %将高频信号处理成间断性
    if mod(t(i),0.25)>0.11&&mod(t(i),0.25)<0.12
    else
        y(i) = 0;
    end
end
sig = x+y; %信号叠加
figure('color','white')
plot(t,sig,'k') %绘制原始信号
%% 2.CEEMDAN分解图
Nstd = 0.2; %Nstd为附加噪声标准差与Y标准差之比
NE = 100;   %NE为对信号的平均次数
MaxIter = 1000;% MaxIter 最大迭代次数
imf = pCEEMDAN(sig,t,Nstd,NE,MaxIter);
% function imf = pCEEMDAN(data,FsOrT,Nstd,NE,MaxIter)
% 画信号CEEMDAN分解图
% 输入：
% data为待分解信号
% FsOrT为采样频率或采样时间向量，如果为采样频率，该变量输入单个值；如果为时间向量，该变量为与y相同长度的一维向量。如果未知采样频率，可设置为1
% Nstd为附加噪声标准差与Y标准差之比
% NE为对信号的平均次数
% MaxIter：最大筛选迭代次数
% 输出：
% imf为经CEEMDAN分解后的各imf分量值
% 例1：（FsOrT为采样频率）
% fs = 100;
% t = 1/fs:1/fs:1;
% data = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDAN(data,fs,0.2,100);
% 例2：（FsOrT为时间向量，需要注意此时FsOrT的长度要与y相同）
% t = 0:0.01:1;
% data = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDAN(data,t,0.2,100);
%% 3.CEEMDAN分解及频谱图
imf = pCEEMDANandFFT(sig,fs,Nstd,NE,MaxIter);% 画信号EEMD分解与各IMF分量频谱对照图
% function imf = pCEEMDANandFFT(y,FsOrT,Nstd,NE,MaxIter)
% 画信号CEEMDAN分解与各IMF分量频谱对照图
% 输入：
% y为待分解信号
% FsOrT为采样频率或采样时间向量，如果为采样频率，该变量输入单个值；如果为时间向量，该变量为与y相同长度的一维向量
% Nstd为附加噪声标准差与Y标准差之比
% NE为对信号的平均次数
% MaxIter：最大筛选迭代次数
% 输出：
% imf为经CEEMDAN分解后的各imf分量值
% 例1：（FsOrT为采样频率）
% fs = 100;
% t = 1/fs:1/fs:1;
% y = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDANandFFT(y,fs,0.2,100);
% 例2：（FsOrT为时间向量，需要注意此时FsOrT的长度要与y相同）
% t = 0:0.01:1;
% y = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDANandFFT(y,t,0.2,100);
