%CEEMDAN�����Ĳ��Խű�����
%% 1.���ɷ����ź�
fs = 400;  %����Ƶ��
t = 0:1/fs:0.75; %ʱ����
x = sin(2*pi*4*t); %��Ƶ�����ź�
y = 0.5*sin(2*pi*120*t); %��Ƶ�����ź�
for i = 1:length(t) %����Ƶ�źŴ���ɼ����
    if mod(t(i),0.25)>0.11&&mod(t(i),0.25)<0.12
    else
        y(i) = 0;
    end
end
sig = x+y; %�źŵ���
figure('color','white')
plot(t,sig,'k') %����ԭʼ�ź�
%% 2.CEEMDAN�ֽ�ͼ
Nstd = 0.2; %NstdΪ����������׼����Y��׼��֮��
NE = 100;   %NEΪ���źŵ�ƽ������
MaxIter = 1000;% MaxIter ����������
imf = pCEEMDAN(sig,t,Nstd,NE,MaxIter);
% function imf = pCEEMDAN(data,FsOrT,Nstd,NE,MaxIter)
% ���ź�CEEMDAN�ֽ�ͼ
% ���룺
% dataΪ���ֽ��ź�
% FsOrTΪ����Ƶ�ʻ����ʱ�����������Ϊ����Ƶ�ʣ��ñ������뵥��ֵ�����Ϊʱ���������ñ���Ϊ��y��ͬ���ȵ�һά���������δ֪����Ƶ�ʣ�������Ϊ1
% NstdΪ����������׼����Y��׼��֮��
% NEΪ���źŵ�ƽ������
% MaxIter�����ɸѡ��������
% �����
% imfΪ��CEEMDAN�ֽ��ĸ�imf����ֵ
% ��1����FsOrTΪ����Ƶ�ʣ�
% fs = 100;
% t = 1/fs:1/fs:1;
% data = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDAN(data,fs,0.2,100);
% ��2����FsOrTΪʱ����������Ҫע���ʱFsOrT�ĳ���Ҫ��y��ͬ��
% t = 0:0.01:1;
% data = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDAN(data,t,0.2,100);
%% 3.CEEMDAN�ֽ⼰Ƶ��ͼ
imf = pCEEMDANandFFT(sig,fs,Nstd,NE,MaxIter);% ���ź�EEMD�ֽ����IMF����Ƶ�׶���ͼ
% function imf = pCEEMDANandFFT(y,FsOrT,Nstd,NE,MaxIter)
% ���ź�CEEMDAN�ֽ����IMF����Ƶ�׶���ͼ
% ���룺
% yΪ���ֽ��ź�
% FsOrTΪ����Ƶ�ʻ����ʱ�����������Ϊ����Ƶ�ʣ��ñ������뵥��ֵ�����Ϊʱ���������ñ���Ϊ��y��ͬ���ȵ�һά����
% NstdΪ����������׼����Y��׼��֮��
% NEΪ���źŵ�ƽ������
% MaxIter�����ɸѡ��������
% �����
% imfΪ��CEEMDAN�ֽ��ĸ�imf����ֵ
% ��1����FsOrTΪ����Ƶ�ʣ�
% fs = 100;
% t = 1/fs:1/fs:1;
% y = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDANandFFT(y,fs,0.2,100);
% ��2����FsOrTΪʱ����������Ҫע���ʱFsOrT�ĳ���Ҫ��y��ͬ��
% t = 0:0.01:1;
% y = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pCEEMDANandFFT(y,t,0.2,100);
