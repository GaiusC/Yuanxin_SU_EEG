function imf = kCEEMDAN(data,Nstd,NE,MaxIter)

% 
% ���룺
% data�����ֽ������
% Nstd��Ϊ����������׼����Y��׼��֮��
% NE�����źŵ�ƽ������
% MaxIter������������
% �����
% imf���ں�ģ̬������ͳһΪn*m��ʽ������nΪģ̬����mΪ���ݵ��������� imf(1,:)��IMF1��imf(emd,:)��Ϊ�в�




%%
[imf]=ceemdan(data,Nstd,NE,MaxIter); %����,MaxIterֱ�����ã����Ը���
end
%%