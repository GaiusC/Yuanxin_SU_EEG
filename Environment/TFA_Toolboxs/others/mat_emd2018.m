function [varargout] = mat_emd2018(x, varargin)
%EMD Empirical Mode Decomposition to 1D signal
%   [IMF, RESIDUAL] = EMD(X) returns intrinsic mode functions (IMF) and
%   residual signal corresponding to the empirical mode decomposition of X,
%   with default settings. X could be a vector or a single data column
%   timetable. When X is a vector, IMF is a 2D matrix where each data
%   column stores an extracted intrinsic mode function, and RESIDUAL is a
%   column vector storing the residual. When X is a single data column
%   timetable, the returned IMF is a timetable and RESIDUAL is a single
%   data column timetable.
%
%   [IMF, RESIDUAL] = EMD(X, 'Name1', Value1, 'Name2', Value2,...)
%   specifies empirical mode decomposition parameters that configure the
%   interpolation method, sifting stop criterion, decomposition stop
%   criterion and boundary methods to be used for empirical mode
%   decomposition. The supported Name-Value Pairs are:
%
%      'SiftRelativeTolerance': One of the sifting stop criteria, also
%                               called Cauchy type convergence criterion,
%                               computed as
%                               ||{c_{i-1}(t)-c_{i}(t)||^2/||c_{i}(t)||^2
%                               at i th Sifting iteration. The default
%                               value is 0.2 [2]. Sifting will stop when
%                               current tolerance is less than
%                               SiftRelativeTolerance.
%
%      'SiftMaxIterations': 	One of the sifting stop criteria, maximum
%                               number of sifting. Default as 100. Sifting
%                               will stop when current iteration is larger
%                               than SiftMaxIterations.
%
%      'MaxNumIMF':           	One of the decomposition stop criteria.
%                               Maximum number of IMF extracted. Default as
%                               10. EMD will stop when number of IMF is
%                               larger than MaxNumIMF.
%
%      'MaxNumExtrema':        	One of the decomposition stop criteria.
%                               Maximum number of extrema in the residual
%                               signal. Default as 1. EMD will stop when
%                               number of extrema is less than
%                               MaxNumExtrema.
%
%      'MaxEnergyRatio':       	One of the decomposition stop criteria.
%                               Signal to residual energy ratio, computed
%                               as 10*log10(||X(t)||^2/||r_k(t)||^2), where
%                               X(t) is the original signal, r_k(t) is the
%                               residual of the k th IMF. The default value
%                               is 20. EMD will stop when energy ratio is
%                               larger than MaxEnergyRatio.
%
%      'Interpolation':        	interpolation method for constructing the
%                               envelope. The method can be one of the
%                               follows: 
%                               'spline' - cubic spline (default)
%                               'pchip'  - piecewise cubic Hermite
%                                          interpolating polynomial, for
%                                          non-smooth signals
%
%      'Display':              	Display outputs for each IMF during
%                               decomposition process (1/0)
%
%   [IMF, RESIDUAL, INFO] = EMD(X, 'Name1', Value1, 'Name2', Value2,...)
%   returns additional information on IMFs and residual for diagnostics
%   purpose.
%
%   EMD(...) with no output arguments plots the original signal, IMFs and
%   residual signal in the current figure.
%
%   EXAMPLE 1: Empirical mode decomposition of non-stationary signals
%   fs = 1000;
%   t = (0:1/fs:4)';
%   x1 = sin(2*pi*50*t) + sin(2*pi*200*t) + 0.1*randn(length(t),1);
%   x2 = sin(2*pi*50/2*t) + sin(2*pi*200/2*t) + sin(2*pi*250*t) + 0.1*randn(length(t),1);
%   x = [x1;x2];
%   [imf, residual, info] = emd(x);
%
%   See also HHT.
%

% Copyright 2017-2018 The MathWorks, Inc.

%#codegen
narginchk(1,inf);
nargoutchk(0,3);

[x,t,td,isTT,opt] = parseAndValidateInputs(x, varargin{:});
[IMF, residual, info] = localEMD(x,t,opt);

if(isTT)
    t = td;
end

if (nargout == 0 && coder.target('MATLAB'))
    signalwavelet.internal.guis.plot.emdPlot(x, IMF, residual, t);
end

if(isTT && coder.target('MATLAB'))
    IMF = array2timetable(IMF,'RowTimes',t);
    residual = array2timetable(residual,'RowTimes',t);
end

if nargout > 0
    varargout{1} = IMF;
end

if nargout > 1
    varargout{2} = residual;
end

if nargout > 2
    varargout{3} = info;
end

end

%--------------------------------------------------------------------------
function [x,t,td,isTT,opt] = parseAndValidateInputs(x,  varargin)
% input type checking
validateattributes(x, {'single','double','timetable'},{'vector'},'emd','X');

% handle timetable and single
isTT = isa(x,'timetable');
isSingle = isa(x,'single');
if(isTT)
    signalwavelet.internal.util.utilValidateattributesTimetable(x, {'sorted','multichannel'},'emd','X');
    [x, t, td] = signalwavelet.internal.util.utilParseTimetable(x);
else
    td = [];
    if(isSingle)
        t = single(1:length(x))';
    else
        t = (1:length(x))';
    end
end

% data integrity checking
validateattributes(x, {'single','double'},{'nonnan','finite','real','nonsparse'},'emd','X');
validateattributes(t, {'single','double'},{'nonnan','finite','real'},'emd','T');

% turn x into column vector
x = x(:);

% parse and validate name-value pairs
if(isempty(varargin))
    opt = signalwavelet.internal.emd.emdOptions();
else
    opt = signalwavelet.internal.emd.emdOptions(varargin{:});
end

validatestring(opt.Interpolation,{'spline', 'pchip'}, 'emd', 'Interpolation');
validateattributes(opt.SiftStopCriterion.SiftMaxIterations,...
    {'numeric'},{'nonnan','finite','scalar','>',0,'integer'}, 'emd', 'SiftMaxIterations');
validateattributes(opt.SiftStopCriterion.SiftRelativeTolerance,...
    {'numeric'},{'nonnan','finite','scalar','>=',0,'<=',1},'emd', 'SiftRelativeTolerance');
validateattributes(opt.DecompStopCriterion.MaxEnergyRatio,...
    {'numeric'},{'nonnan','finite','scalar'}, 'emd', 'MaxEnergyRatio');
validateattributes(opt.DecompStopCriterion.MaxNumExtrema,...
    {'numeric'},{'nonnan','finite','scalar','>=',0,'integer'},'emd','MaxNumExtrema');
validateattributes(opt.DecompStopCriterion.MaxNumIMF,...
    {'numeric'},{'nonnan','finite','scalar','>',0,'integer'},'emd', 'MaxNumIMF');
end

%--------------------------------------------------------------------------
function [IMFs, rsig, info] = localEMD(x, t, opt)
isInMATLAB = coder.target('MATLAB');
isSingle = isa(x,'single');

% get name-value pairs
Interpolation = opt.Interpolation;
MaxEnergyRatio = opt.DecompStopCriterion.MaxEnergyRatio;
MaxNumExtrema = opt.DecompStopCriterion.MaxNumExtrema;
MaxNumIMF = opt.DecompStopCriterion.MaxNumIMF;
SiftMaxIterations = opt.SiftStopCriterion.SiftMaxIterations;
SiftRelativeTolerance = opt.SiftStopCriterion.SiftRelativeTolerance;
Display = opt.Display;

% initialization
rsig = x;
N = length(x);

if(isSingle)
    ArrayType = 'single';
else
    ArrayType = 'double';
end

IMFs = zeros(N, MaxNumIMF, ArrayType);
info.NumIMF = zeros(MaxNumIMF, 1, ArrayType);
info.NumExtrema = zeros(MaxNumIMF, 1, ArrayType);
info.NumZerocrossing = zeros(MaxNumIMF, 1, ArrayType);
info.NumSifting = zeros(MaxNumIMF, 1, ArrayType);
info.MeanEnvelopeEnergy = zeros(MaxNumIMF, 1, ArrayType);
info.RelativeTolerance = zeros(MaxNumIMF, 1, ArrayType);

% preallocate memory
rsigPrev = zeros(N, 1, ArrayType);
mVal = zeros(N, 1, ArrayType);
upperEnvelope = zeros(N, 1, ArrayType);
lowerEnvelope = zeros(N, 1, ArrayType);

% Define intermediate print formats
if(isInMATLAB && opt.Display)
    fprintf('Current IMF  |  #Sift Iter  |  Relative Tol  |  Stop Criterion Hit  \n');
    formatstr = '  %5.0f      |    %5.0f     | %12.5g   |  %s\n';
end

% use different functions under different environment
if(isInMATLAB)
    if(~isSingle)
        localFindExtramaIdx = @(x) signalwavelet.internal.emd.cg_utilFindExtremaIdxmex_double(x);
    else
        localFindExtramaIdx = @(x) signalwavelet.internal.emd.cg_utilFindExtremaIdxmex_single(x);
    end
else
    localFindExtramaIdx = @(x) signalwavelet.internal.emd.utilFindExtremaIdx(x);
end

% extract IMFs
i = 0;
while(i<MaxNumIMF)
    % convergence checking
    [peaksIdx, bottomsIdx] = localFindExtramaIdx(rsig);
    numResidExtrema = length(peaksIdx) + length(bottomsIdx);
    energyRatio = 10*log10(norm(x,2)/norm(rsig,2));
    
    if(energyRatio>MaxEnergyRatio || numResidExtrema<MaxNumExtrema)
        break
    end
    
    % SIFTING process initialization
    rsigL = rsig;
    rtol = ones(1, ArrayType);
    k = 0;
    SiftStopCriterionHit = 'SiftMaxIteration';
    
    % Sifting process
    while (k<SiftMaxIterations)
        % check convergence
        if(rtol<SiftRelativeTolerance)
            SiftStopCriterionHit = 'SiftMaxRelativeTolerance';
            break;
        end
        
        % store previous residual
        rsigPrev(1:N) = rsigL;
        
        % finding peaks
        [peaksIdx, bottomsIdx] = localFindExtramaIdx(rsigL);
        
        if((length(peaksIdx) + length(bottomsIdx))>0)
            % compute upper and lower envelope using extremas
            [uLoc, uVal, bLoc, bVal] = computeSupport(t, rsigL, peaksIdx, bottomsIdx);
            upperEnvelope(:) = interp1(uLoc, uVal, t, Interpolation);
            lowerEnvelope(:) = interp1(bLoc, bVal, t, Interpolation);
            
            % subtract mean envelope from residual
            mVal(1:N) = (upperEnvelope + lowerEnvelope)/2;
        else
            mVal(1:N) = 0;
        end
        
        rsigL = rsigL - mVal;
        
        % residual tolerance
        rtol = (norm(rsigPrev-rsigL,2)/norm(rsigPrev,2))^2;
        k = k + 1;
    end
    
    if(isInMATLAB && Display)
        fprintf(formatstr, i+1, k, rtol, SiftStopCriterionHit);
    end
    
    % record information
    [peaksIdx, bottomsIdx] = localFindExtramaIdx(rsigL);
    numZerocrossing = sum(diff(sign(rsigL))~=0);
    info.NumIMF(i+1) = i+1;
    info.NumExtrema(i+1) = length(peaksIdx) + length(bottomsIdx);
    info.NumZerocrossing(i+1) = numZerocrossing;
    info.MeanEnvelopeEnergy(i+1) = mean(mVal.^2);
    info.NumSifting(i+1) = k;
    info.RelativeTolerance(i+1) = rtol;
    
    % after sifting
    if(k==0)
        break;
    end
    
    % extract new IMF and subtract the IMF from residual signal
    IMFs(:,i+1) = rsigL;
    rsig = rsig - IMFs(:,i+1);
    i = i + 1;
end

% remove extra portion
IMFs = IMFs(:,1:i);
info.NumIMF = info.NumIMF(1:i);
info.NumExtrema = info.NumExtrema(1:i);
info.NumZerocrossing = info.NumZerocrossing(1:i);
info.NumSifting = info.NumSifting(1:i);
info.MeanEnvelopeEnergy = info.MeanEnvelopeEnergy(1:i);
info.RelativeTolerance = info.RelativeTolerance(1:i);
end

%--------------------------------------------------------------------------
function [uLoc, uVal, bLoc, bVal] = computeSupport(t, rsigL, pksIdx, btmsIdx)
% compute support for upper and lower envelope given input signal rsigL
N = length(t);
if(isempty(pksIdx))
    pksIdx = [1; N];
end

if(isempty(btmsIdx))
    btmsIdx = [1; N];
end

pksLoc = t(pksIdx);
btmsLoc = t(btmsIdx);

% compute envelop for wave method
% extended waves on the left
[lpksLoc, lpksVal, lbtmLoc, lbtmVal] = signalwavelet.internal.emd.emdWaveExtension(t(1), rsigL(1),...
    pksLoc(1), rsigL(pksIdx(1)),...
    btmsLoc(1), rsigL(btmsIdx(1)),...
    -1);

% extended waves on the right
[rpksLoc, rpksVal, rbtmLoc, rbtmVal] = signalwavelet.internal.emd.emdWaveExtension(t(end), rsigL(end),...
    pksLoc(end), rsigL(pksIdx(end)),...
    btmsLoc(end), rsigL(btmsIdx(end)),...
    1);

% append extended wave to extrema
uLoc = [lpksLoc;pksLoc;rpksLoc];
uVal = [lpksVal;rsigL(pksIdx);rpksVal];
bLoc = [lbtmLoc;btmsLoc;rbtmLoc];
bVal = [lbtmVal;rsigL(btmsIdx);rbtmVal];
end

