function imf = pCEEMDANandFFT(y,FsOrT,Nstd,NE,MaxIter)
% 画信号CEEMDAN分解与各IMF分量频谱对照图
% 输入：
% y为待分解信号
% FsOrT为采样频率或采样时间向量，如果为采样频率，该变量输入单个值；如果为时间向量，该变量为与y相同长度的一维向量
% Nstd为附加噪声标准差与Y标准差之比
% NE为对信号的平均次数
% MaxIter：最大迭代次数
% 输出：
% imf为经CEEMDAN分解后的各imf分量值
% 例1：（FsOrT为采样频率）
% fs = 100;
% t = 1/fs:1/fs:1;
% y = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDANandFFT(y,fs,0.2,100,1000);
% 例2：（FsOrT为时间向量，需要注意此时FsOrT的长度要与y相同）
% t = 0:0.01:1;
% y = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDANandFFT(y,t,0.2,100,1000);

% 注意：在使用该代码之前，请务必安装工具箱：http://www.khscience.cn/docs/index.php/2020/04/09/1/

%  Copyright (c) 2020 Mr.括号 All rights reserved.
%  本代码为淘宝买家专用，不开源，请勿公开分享~
%%

if length(FsOrT) == 1
    t = 1/FsOrT:1/FsOrT:length(y)/FsOrT;
    Fs = FsOrT;
else
    t = FsOrT;
    Fs = 1/(t(2)-t(1));
end
imf = kCEEMDAN(y,Nstd,NE,MaxIter);
figure('Name','CEEMDAN分解与各IMF分量频谱对照图','Color','white');
subplot(size(imf,1)+1,2,1);
plot(t,y,'k');grid on;
ylabel('原始数据');
title('CEEMDAN分解');
set(gca,'XTick',[]);
subplot(size(imf,1)+1,2,2);
pFFT(y,Fs);grid on;
title('对应频谱');
set(gca,'XTick',[]);
for i = 2:size(imf,1)+1
    subplot(size(imf,1)+1,2,i*2-1);
    plot(t,imf(i-1,:),'k');
    ylabel(['IMF',num2str(i-1)]);
    if (i ~= size(imf,1)+1)
        set(gca,'XTick',[]);
    end
    if (i == size(imf,1)+1)
        ylabel(['res']);
        xlabel('time/s');
    end
    grid on;
    subplot(size(imf,1)+1,2,i*2);
    pFFT(imf(i-1,:),Fs);
    if (i ~= size(imf,1)+1)
        set(gca,'XTick',[]);
    end
    if (i == size(imf,1)+1)
        xlabel('frequency/Hz');
    end
    grid on;
end
end
