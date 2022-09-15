function imf = pCEEMDAN(data,FsOrT,Nstd,NE,MaxIter)
% 画信号CEEMDAN分解图
% 输入：
% data为待分解信号
% FsOrT为采样频率或采样时间向量，如果为采样频率，该变量输入单个值；如果为时间向量，该变量为与y相同长度的一维向量。如果未知采样频率，可设置为1
% Nstd为附加噪声标准差与Y标准差之比
% NE为对信号的平均次数
% MaxIter：最大迭代次数
% 输出：
% imf为经CEEMDAN分解后的各imf分量值
% 例1：（FsOrT为采样频率）
% fs = 100;
% t = 1/fs:1/fs:1;
% data = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDAN(data,fs,0.2,100,1000);
% 例2：（FsOrT为时间向量，需要注意此时FsOrT的长度要与y相同）
% t = 0:0.01:1;
% data = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDAN(data,t,0.2,100,1000);

%
if length(FsOrT) == 1  %如果输入的为频率值
    t = 1/FsOrT:1/FsOrT:length(data)/FsOrT;
else
    t = FsOrT;         %如果输入的为时间向量
end
imf=kCEEMDAN(data,Nstd,NE,MaxIter);
rows = size(imf,1);    %获取分量数目

figure('Name','CEEMDAN','Color','white');
subplot(rows+1,1,1);
plot(t,data);grid on;
xlim([t(1) t(end)]);
ylabel('Raw data');
title('CEEMDAN');

for i = 1:size(imf,1)
    subplot(rows+1,1,i+1);
    plot(t,imf(i,:));
    xlim([t(1) t(end)]);
    ylabel(['IMF',num2str(i)]);
    if (i == size(imf,1))
        ylabel(['res']);
        xlabel('time');
    end
    grid on;
end
end

