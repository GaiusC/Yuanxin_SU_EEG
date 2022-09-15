function [A,W,icasig,x] = ICA_OD(x,input_row,wlength,alfa)
% for n = 1:input_row
    [icasig, A, W] = fastica(x,'numOfIC',input_row,'displayMode','off','firstEig',1,'lastEig',1,'g','gauss');
%     'firstEig',1,'lastEig',input_row
    [row,col] = size(icasig);
    alfa_array = [1.38,1.54,1.65,1.73,1.8,1.86,1.92,1.96];
    
    
    for k = 1:row
%         TF = isoutlier(icasig(k,:),'movmean',10);
%         icasig(k,TF) = 0;
        icmean = mean(icasig(k,:));
        icstd = std(icasig(k,:));
        dmax = icstd*alfa;
        for l = 1:wlength:col-wlength
            for i = 1:wlength
                di = abs(icasig(k,l+i-1) - icmean);
                if di <= dmax
                    icasig(k,l+i-1) =0;% median(icasig(k,:));
                end
            end
        end
    end
    x=A*icasig;
end
% end