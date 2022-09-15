function imf = pCEEMDANandFFT(y,FsOrT,Nstd,NE,MaxIter)
% ���ź�CEEMDAN�ֽ����IMF����Ƶ�׶���ͼ
% ���룺
% yΪ���ֽ��ź�
% FsOrTΪ����Ƶ�ʻ����ʱ�����������Ϊ����Ƶ�ʣ��ñ������뵥��ֵ�����Ϊʱ���������ñ���Ϊ��y��ͬ���ȵ�һά����
% NstdΪ����������׼����Y��׼��֮��
% NEΪ���źŵ�ƽ������
% MaxIter������������
% �����
% imfΪ��CEEMDAN�ֽ��ĸ�imf����ֵ
% ��1����FsOrTΪ����Ƶ�ʣ�
% fs = 100;
% t = 1/fs:1/fs:1;
% y = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDANandFFT(y,fs,0.2,100,1000);
% ��2����FsOrTΪʱ����������Ҫע���ʱFsOrT�ĳ���Ҫ��y��ͬ��
% t = 0:0.01:1;
% y = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDANandFFT(y,t,0.2,100,1000);

% ע�⣺��ʹ�øô���֮ǰ������ذ�װ�����䣺http://www.khscience.cn/docs/index.php/2020/04/09/1/

%  Copyright (c) 2020 Mr.���� All rights reserved.
%  ������Ϊ�Ա����ר�ã�����Դ�����𹫿�����~
%%

if length(FsOrT) == 1
    t = 1/FsOrT:1/FsOrT:length(y)/FsOrT;
    Fs = FsOrT;
else
    t = FsOrT;
    Fs = 1/(t(2)-t(1));
end
imf = kCEEMDAN(y,Nstd,NE,MaxIter);
figure('Name','CEEMDAN�ֽ����IMF����Ƶ�׶���ͼ','Color','white');
subplot(size(imf,1)+1,2,1);
plot(t,y,'k');grid on;
ylabel('ԭʼ����');
title('CEEMDAN�ֽ�');
set(gca,'XTick',[]);
subplot(size(imf,1)+1,2,2);
pFFT(y,Fs);grid on;
title('��ӦƵ��');
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
