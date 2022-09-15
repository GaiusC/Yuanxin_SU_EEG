function imf = pCEEMDAN(data,FsOrT,Nstd,NE,MaxIter)
% ���ź�CEEMDAN�ֽ�ͼ
% ���룺
% dataΪ���ֽ��ź�
% FsOrTΪ����Ƶ�ʻ����ʱ�����������Ϊ����Ƶ�ʣ��ñ������뵥��ֵ�����Ϊʱ���������ñ���Ϊ��y��ͬ���ȵ�һά���������δ֪����Ƶ�ʣ�������Ϊ1
% NstdΪ����������׼����Y��׼��֮��
% NEΪ���źŵ�ƽ������
% MaxIter������������
% �����
% imfΪ��CEEMDAN�ֽ��ĸ�imf����ֵ
% ��1����FsOrTΪ����Ƶ�ʣ�
% fs = 100;
% t = 1/fs:1/fs:1;
% data = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDAN(data,fs,0.2,100,1000);
% ��2����FsOrTΪʱ����������Ҫע���ʱFsOrT�ĳ���Ҫ��y��ͬ��
% t = 0:0.01:1;
% data = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDAN(data,t,0.2,100,1000);

%
if length(FsOrT) == 1  %��������ΪƵ��ֵ
    t = 1/FsOrT:1/FsOrT:length(data)/FsOrT;
else
    t = FsOrT;         %��������Ϊʱ������
end
imf=kCEEMDAN(data,Nstd,NE,MaxIter);
rows = size(imf,1);    %��ȡ������Ŀ

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

