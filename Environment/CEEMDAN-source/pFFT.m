function pFFT(y,Fs)
%  �ú������������룬���ź�ֵ�����Ƶ��
%  �ú���ֻ��һ��ͼ������plot����δʹ��figure
%  Copyright (c) 2020 Mr.���� All rights reserved.
%  ������Ϊ�Ա����ר�ã�����Դ�����𹫿�����~
%% 
t_s = 1/Fs; %��������
t_start = 0;
t_end = (length(y)-1)*t_s;
t = 0 : t_s : t_end;
Druation = t_end -t_start;  %�������ʱ��
Sampling_points = Druation/t_s +1;  %��������
nfft = Sampling_points;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ƶ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_f = fft(y,nfft); %����Ҷ�任
f_s = 1/t_s; %����Ƶ��
f_x = 0:f_s/(Sampling_points -1):f_s;  %ע������ͺ�����Ƶ�ʶ�Ӧ���ˣ�Ƶ�ʷֱ��ʾ���f_s/(Sampling_points -1)
plot(f_x(1:round(length(f_x)/2)),2/Sampling_points*abs(y_f(1:round(length(f_x)/2))),'k');
