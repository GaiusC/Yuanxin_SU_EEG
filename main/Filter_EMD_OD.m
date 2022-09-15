%% same table and figure show the different measurements of different algorithm and 

clc,clear,close all
warning off
%%
% Dataset = readmatrix('C:\Desktop\EEG-EyeBlinks/EEG-IO/S01_data.csv');
%label = readmatrix('/Users/chenyiya/Desktop/EEG-EyeBlinks/EEG-IO/S01_labels.csv');
% Original_Data = Dataset(1:length(time),3)'/1e3;

%% Semi_data
Pure_Data = load('Pure_Data.mat');
Contaminated_Data = load('Contaminated_Data.mat');
Original_Data = Contaminated_Data.sim1_con(2,1:3000);
Pure_Data_C1 = Pure_Data.sim1_resampled(2,1:3000);

%% parameter initialization
fs = 200;
time = 0:1/fs:(length(Original_Data)-1)/fs;
Nstd = 0.2;%Noise standard deviation
NE = 8;
MaxIter=10;%maximum times of iteration
% lev=7;%Level of wavelet decomposition
alfa = 1.92;
wlength = 9;


%% Data Decomposition
LPfilter = LowPass_4Hz_Hamming();
% Original_Data_LP = conv(Original_Data,LPfilter.Numerator,'same');
% Original_Data_HP = Original_Data-Original_Data_LP;
CEEMDAN = pCEEMDAN(Original_Data,fs,Nstd,NE,MaxIter);
EMD = emd(Original_Data,'Display',1)';
EEMD = eemd(Original_Data,Nstd,NE)';
EEMD(1,:)=[];
% CEEMDAN(end-5,:)=[];
[CEEMDAN_row,~]=size(CEEMDAN);
[EMD_row,~]=size(EMD);
[EEMD_row,~]=size(EEMD);
%% CEEMDAN-FILTER-ICA-OA

parfor n = 1:CEEMDAN_row
    CEEMDAN_LP(n,:) = conv(CEEMDAN(n,:),LPfilter.Numerator,'same');
end
parfor n = 1:EMD_row
    EMD_LP(n,:) = conv(EMD(n,:),LPfilter.Numerator,'same');
end
parfor n = 1:EEMD_row
    EEMD_LP(n,:) = conv(EEMD(n,:),LPfilter.Numerator,'same');
end
EMD_HP = EMD-EMD_LP;
EEMD_HP = EEMD - EEMD_LP;
CEEMDAN_HP = CEEMDAN - CEEMDAN_LP;
[EMD_row,EMD_col]=size(EMD_LP);
[EEMD_row,EEMD_col]=size(EEMD_LP);
[CEEMDAN_row,CEEMDAN_col]=size(CEEMDAN_LP);


[EMD_Th,EMD_Non_Th,EMD_count_Th,EMD_count_Non_Th,EMD_PSD,EMD_freq,EMD_Entropy,EMD_I_Entropy,EMD_H] = PSD_Thesholding(EMD_LP,EMD_row,'EMD');
[EEMD_Th,EEMD_Non_Th,EEMD_count_Th,EEMD_count_Non_Th,EEMD_PSD,EEMD_freq,EEMD_Entropy,EEMD_I_Entropy,EEMD_H] = PSD_Thesholding(EEMD_LP,EEMD_row,'EEMD');
[CEEMDAN_Th,CEEMDAN_Non_Th,CEEMDAN_count_Th,CEEMDAN_count_Non_Th,CEEMDAN_PSD,CEEMDAN_freq,CEEMDAN_Entropy,CEEMDAN_I_Entropy,CEEMDAN_H] = PSD_Thesholding(CEEMDAN_LP,CEEMDAN_row,'CEEMDAN');


[EMD_A,W,EMD_icasig,EMD_Th] = ICA_OD(EMD_Th,EMD_count_Th,wlength,alfa);
[EEMD_A,W,EEMD_icasig,EEMD_Th] = ICA_OD(EEMD_Th,EEMD_count_Th,wlength,alfa);
[CEEMDAN_A,W,CEEMDAN_icasig,CEEMDAN_Th] = ICA_OD(CEEMDAN_Th,CEEMDAN_count_Th,wlength,alfa);


EMD_FIOD_Result = sum(EMD_Th)+sum(EMD_HP);
EEMD_FIOD_Result = sum(EEMD_Th)+sum(EEMD_HP);
CEEMDAN_FIOD_Result = sum(CEEMDAN_Th)+sum(CEEMDAN_HP);

EMD_Rem_part =  Original_Data-EMD_FIOD_Result;
EEMD_Rem_part = Original_Data-EEMD_FIOD_Result;
CEEMDAN_Rem_part = Original_Data-CEEMDAN_FIOD_Result;


Ori_Noise = Original_Data-Pure_Data_C1;

%% PSD Calculate
[EMD_PSD_Rem,EMD_freq_Rem] = periodogram(EMD_Rem_part,[],[],fs);
[EMD_PSD_result,EMD_freq_result] = periodogram(EMD_FIOD_Result,[],[],fs);

[EEMD_PSD_Rem,EEMD_freq_Rem] = periodogram(EEMD_Rem_part,[],[],fs);
[EEMD_PSD_result,EEMD_freq_result] = periodogram(EEMD_FIOD_Result,[],[],fs);

[CEEMDAN_PSD_Rem,CEEMDAN_freq_Rem] = periodogram(CEEMDAN_Rem_part,[],[],fs);
[CEEMDAN_PSD_result,CEEMDAN_freq_result]=periodogram(CEEMDAN_FIOD_Result,[],[],fs);

[PSD_OriN,freq_OriN]=periodogram(Ori_Noise,[],[],fs);
[PSD_Pure,freq_Pure] = periodogram(Pure_Data_C1,[],[],fs);

%% DISPLAY
EMD_RRMSE = rms(EMD_FIOD_Result-Pure_Data_C1,1)/rms(Pure_Data_C1,1);
EEMD_RRMSE = rms(EEMD_FIOD_Result-Pure_Data_C1,1)/rms(Pure_Data_C1,1);
CEEMDAN_RRMSE = rms(CEEMDAN_FIOD_Result-Pure_Data_C1,1)/rms(Pure_Data_C1,1);

EMD_RRMSE_Noise = rms(Ori_Noise-EMD_Rem_part,1)/rms(Ori_Noise,1);
EEMD_RRMSE_Noise = rms(Ori_Noise-EEMD_Rem_part,1)/rms(Ori_Noise,1);
CEEMDAN_RRMSE_Noise = rms(Ori_Noise-CEEMDAN_Rem_part,1)/rms(Ori_Noise,1);

EMD_RRMSE_PSD_Noise = rms(PSD_OriN-EMD_PSD_Rem,1)/rms(PSD_OriN,1);
EEMD_RRMSE_PSD_Noise = rms(PSD_OriN-EEMD_PSD_Rem,1)/rms(PSD_OriN,1);
CEEMDAN_RRMSE_PSD_Noise = rms(PSD_OriN-CEEMDAN_PSD_Rem,1)/rms(PSD_OriN,1);

% EMD_SNRv = sqrt(sum((EMD_FIOD_Result/1e6).^2)/sum((EMD_Rem_part/1e6).^2));
% EMD_SNR = 20*log10(EMD_SNRv);
% 
% EEMD_SNRv = sqrt(sum((EEMD_FIOD_Result/1e6).^2)/sum((EEMD_Rem_part/1e6).^2));
% EEMD_SNR = 20*log10(EEMD_SNRv);
% 
% CEEMDAN_SNRv = sqrt(sum((CEEMDAN_FIOD_Result/1e6).^2)/sum((CEEMDAN_Rem_part/1e6).^2));
% CEEMDAN_SNR = 20*log10(CEEMDAN_SNRv);

name = {'EMD';'EEMD';'CEEMDAN'};
RRMSE = {string(EMD_RRMSE);string(EEMD_RRMSE);string(CEEMDAN_RRMSE)};
Noise_RRMSE = {string(EMD_RRMSE_Noise);string(EEMD_RRMSE_Noise);string(CEEMDAN_RRMSE_Noise)};
Noise_PSD_RRMSE = {string(EMD_RRMSE_PSD_Noise);string(EEMD_RRMSE_PSD_Noise);string(CEEMDAN_RRMSE_PSD_Noise)};
colname = {'Algorithm';'RRMSE';'RRMSE of Noise';'RRMSE of PSD of Noise'};
Result_Table = table(name,RRMSE,Noise_RRMSE,Noise_PSD_RRMSE,'VariableNames',colname);
disp(Result_Table)
%% Display of result of Decomposition
figure()
for n = 1:EMD_row
    subplot(EMD_row,1,n)
    plot(time,EMD(n,:))
end

figure()
for n = 1:EEMD_row
    subplot(EEMD_row,1,n)
    plot(time,EEMD(n,:))
end
%% Result Compare with Original signal

figure()
title('Result Compare with Original Signal')
subplot(311)
plot(time,Original_Data,time,EMD_FIOD_Result,'r--')
subtitle('EMD','FontSize',20)
legend('Original Signal','Filtered Signal','fontsize',20)
xlabel('time(s)')
ylabel('Amplitude(uV)')
subplot(312)
plot(time,Original_Data,time,EEMD_FIOD_Result,'r--')
subtitle('EEMD','FontSize',20)
legend('Original Signal','Filtered Signal','fontsize',20)
xlabel('time(s)')
ylabel('Amplitude(uV)')
subplot(313)
plot(time,Original_Data,time,CEEMDAN_FIOD_Result,'r--')
subtitle('CEEMDAN','FontSize',20)
legend('Original Signal','Filtered Signal','fontsize',20)
xlabel('time(s)')
ylabel('Amplitude(uV)')

%% Result Compare with Pure signal

figure()
title('Result Compare with Original Signal')
subplot(311)
plot(time,Pure_Data_C1,time,EMD_FIOD_Result,'r--')
subtitle('EMD')
legend('Pure Signal','Filtered Signal by EMD')
xlabel('time(s)')
ylabel('Amplitude(uV)')
subplot(312)
plot(time,Pure_Data_C1,time,EEMD_FIOD_Result,'r--')
subtitle('EEMD')
legend('Pure Signal','Filtered Signal by EEMD')
xlabel('time(s)')
ylabel('Amplitude(uV)')
subplot(313)
plot(time,Pure_Data_C1,time,CEEMDAN_FIOD_Result,'r--')
subtitle('CEEMDAN')
legend('Pure Signal','Filtered Signal')
xlabel('time(s)')
ylabel('Amplitude(uV)')
%% Removal part compare with Original noise
figure()
title('Removal part and Original Noise')
subplot(311)
plot(time,Ori_Noise,'r',time,EMD_Rem_part,'b')
subtitle('EMD','FontSize',20)
legend('Original Noise','Removal part','fontsize',15)
subplot(312)
plot(time,Ori_Noise,'r',time,EEMD_Rem_part,'b')
subtitle('EEMD','FontSize',20)
legend('Original Noise','Removal part','fontsize',15)
subplot(313)
plot(time,Ori_Noise,'r',time,CEEMDAN_Rem_part,'b')
subtitle('CEEMDAN','FontSize',20)
legend('Original Noise','Removal part','fontsize',15)


%% PSD show

figure()
subplot(311)
plot(freq_Pure,PSD_Pure,EMD_freq_result,EMD_PSD_result,'r--','LineWidth',2)
legend('PSD of Pure Signal','PSD of Filtered Signal','FontSize',15,'location','north')
subtitle('EMD','FontSize',20)
xlim([0 10])
xlabel('Frequency(Hz)')
ylabel('Amplitude')
subplot(312)
plot(freq_Pure,PSD_Pure,EEMD_freq_result,EEMD_PSD_result,'r--','LineWidth',2)
legend('PSD of Pure Signal','PSD of Filtered Signal','FontSize',15,'location','north')
subtitle('EEMD','FontSize',20)
xlim([0 10])
xlabel('Frequency(Hz)')
ylabel('Amplitude')
subplot(313)
plot(freq_Pure,PSD_Pure,CEEMDAN_freq_result,CEEMDAN_PSD_result,'r--','LineWidth',2)
legend('PSD of Pure Signal','PSD of Filtered Signal','FontSize',15,'location','north')
subtitle('CEEMDAN','FontSize',20)
xlim([0 10])
xlabel('Frequency(Hz)')
ylabel('Amplitude')

%% PSD of Noise
figure()
subplot(311)
plot(freq_OriN,PSD_OriN,EMD_freq_Rem,EMD_PSD_Rem,'r--','LineWidth',2)
legend('PSD of Original Noise','PSD of Remove part','FontSize',15)
subtitle('EMD','FontSize',20)
xlim([0 8])
xlabel('Frequency(Hz)')
ylabel('Amplitude')
subplot(312)
plot(freq_OriN,PSD_OriN,EEMD_freq_Rem,EEMD_PSD_Rem,'r--','LineWidth',2)
legend('PSD of Original Noise','PSD of Remove part','FontSize',15)
subtitle('EEMD','FontSize',20)
xlim([0 8])
xlabel('Frequency(Hz)')
ylabel('Amplitude')
subplot(313)
plot(freq_OriN,PSD_OriN,CEEMDAN_freq_Rem,CEEMDAN_PSD_Rem,'r--','LineWidth',2)
legend('PSD of Original Noise','PSD of Remove part','FontSize',15)
subtitle('CEEMDAN','FontSize',20)
xlim([0 8])
xlabel('Frequency(Hz)')
ylabel('Amplitude')

figure()
subplot(221)
histogram(Original_Data)
subtitle('PDF of Original Signal')
gtext(strcat('Skewness=',string(skewness(Original_Data))),'Fontsize',20)
gtext(strcat('Kurtosis=',string(kurtosis(Original_Data))),'Fontsize',20)
subplot(222)
histogram(EMD_FIOD_Result)
subtitle('PDF of Filtered Signal by EMD')
gtext(strcat('Skewness=',string(skewness(EMD_FIOD_Result))),'Fontsize',20)
gtext(strcat('Kurtosis=',string(kurtosis(EMD_FIOD_Result))),'Fontsize',20)
subplot(223)
histogram(EEMD_FIOD_Result)
subtitle('PDF of Filtered Signal by EEMD')
gtext(strcat('Skewness=',string(skewness(EEMD_FIOD_Result))),'Fontsize',20)
gtext(strcat('Kurtosis=',string(kurtosis(EEMD_FIOD_Result))),'Fontsize',20)
subplot(224)
histogram(CEEMDAN_FIOD_Result)
subtitle('PDF of Filtered Signal by CEEMDAN')
gtext(strcat('Skewness=',string(skewness(CEEMDAN_FIOD_Result))),'Fontsize',20)
gtext(strcat('Kurtosis=',string(kurtosis(CEEMDAN_FIOD_Result))),'Fontsize',20)

figure()
histogram(Pure_Data_C1)
subtitle('PDF of Original Signal')
gtext(strcat('Skewness=',string(skewness(Pure_Data_C1))),'Fontsize',40)
gtext(strcat('Kurtosis=',string(kurtosis(Pure_Data_C1))),'Fontsize',40)


figure()
plot(time,Original_Data,'LineWidth',2)
xlabel('time(s)','FontSize',40)
ylabel('Amplitude(uV)','FontSize',40)
legend('Comtaminated EEG','fontsize',40)
figure()
plot(time,Pure_Data_C1)
xlabel('time(s)','FontSize',40)
ylabel('Amplitude(uV)','FontSize',40)
legend('Eyeblinking-free EEG','fontsize',40)
figure()
plot(time,Ori_Noise)
xlabel('time(s)','FontSize',40)
ylabel('Amplitude(uV)','FontSize',40)
legend('Ocular Artifacts','fontsize',40)

