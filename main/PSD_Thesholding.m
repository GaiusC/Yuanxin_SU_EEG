function [C_Th,C_Non_Th,count_Th,count_Non_Th,x_PSD,freq,PSD_Entropy,I_Entropy,H] = PSD_Thesholding(x,nrow,method)
C_Non_Th = zeros(1,length(x));
parfor n = 1:nrow
    x(n,:) = (x(n,:)-mean(x(n,:)))./(max(abs(x(n,:))));
    [x_PSD(n,:),freq(n,:)] =periodogram(x(n,:));
    %       [x_PSD(n,:),freq(n,:)] = pwelch(x(n,:),10,[],2048,fs)
end
count_Th=1;
count_Non_Th=1;
for n = 1:nrow
    Energy(n) = sum(x_PSD(n,:).^2);
    Energy_PDF(n,:) = x_PSD(n,:)/Energy(n);
    PSD_Entropy(n) = -sum(Energy_PDF(n,:).*log(Energy_PDF(n,:)));
    I_Entropy(n) = entropy(x(n,:));
    figure(2)
    p = histogram(x(n,:));
    p = p.Values./numel(x(n,:));
    H(n) = -sum(p.*log(p));
    close(2)
end


% s_PSD_Entropy = sort(PSD_Entropy);
% s_I_Entropy = sort(I_Entropy);
if strcmp(method,'EMD')
    PSD_Thresholding = 0.003;
    I_Thresholding = 2.7;%(max(CEEMDAN_Entropy)-min(CEEMDAN_Entropy))/2; PSD=0.005
    %% OD
%     PSD_Thresholding = 2.3;
%     I_Thresholding = 3;
end
if strcmp(method,'EEMD')
    PSD_Thresholding = 0.003;
    I_Thresholding = 2.7;%(max(CEEMDAN_Entropy)-min(CEEMDAN_Entropy))/2; PSD=0.005
    %% OD
%     PSD_Thresholding = 0;
%     I_Thresholding = 3;
end
if strcmp(method,'CEEMDAN')
    PSD_Thresholding = 0.004;
    I_Thresholding = 2.8;%(max(CEEMDAN_Entropy)-min(CEEMDAN_Entropy))/2; PSD=0.005
    %% OD
%     PSD_Thresholding = 2.3;
%     I_Thresholding = 3;
end
% PSD_Entropy = sort(PSD_Entropy);
% Thresholding = PSD_Entropy(6);
for n = 1:nrow-1
    if  ((H(n)) < I_Thresholding) && (PSD_Entropy(n) > PSD_Thresholding)
        C_Th(count_Th,:) = x(n,:);
        count_Th = count_Th+1;
    else
        C_Non_Th(count_Non_Th,:) = x(n,:);
        count_Non_Th= count_Non_Th+1;
        n;
    end
end

% C_Th(count_Th:count_Th+1,:) = x(end-1:end,:);
count_Th = count_Th-1;
count_Non_Th = count_Non_Th-1;
if count_Th == 0
    C_Th = 0;
end
