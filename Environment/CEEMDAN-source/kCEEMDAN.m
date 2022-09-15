function imf = kCEEMDAN(data,Nstd,NE,MaxIter)

% 
% 输入：
% data：待分解的数据
% Nstd：为附加噪声标准差与Y标准差之比
% NE：对信号的平均次数
% MaxIter：最大迭代次数
% 输出：
% imf：内涵模态分量，统一为n*m格式，其中n为模态数，m为数据点数。例如 imf(1,:)即IMF1，imf(emd,:)即为残差




%%
[imf]=ceemdan(data,Nstd,NE,MaxIter); %调用,MaxIter直接设置，可以更改
end
%%