function [varargout] = hht_fix(IMF, varargin)
% 修正过后的hht函数
% 修改了第106行的问题
%HHT Hilbert Spectrum of a signal using Hilbert-Huang Transform
%   P = HHT(IMF) returns the Hilbert Spectrum of the signal specified by
%   IMF (intrinsic mode function) in P. Each column of IMF will be treated
%   as an intrinsic mode function. When IMF is a matrix, normalized
%   frequency will be used as sampling frequency. When IMF is a column or
%   row vector, it is treated as one intrinsic mode function. When IMF is a
%   uniformly-sampled timetable, its sampling frequency Fs will be inferred
%   from the timetable. The output P is a matrix storing the Hilbert
%   Spectrum.
%
%   P = HHT(IMF, fs) specifies the sampling frequency Fs when IMF is a
%   matrix.
%
%   [P,F,T] = HHT(...) returns frequency vector F and time vector T in
%   addition to P.
%
%   [P,F,T,IMFINSF,IMFINSE] = HHT(...) also returns instantaneous
%   frequencies of each IMF in IMFINSF and instantaneous energy of each IMF
%   in IMFINSE. IMFINSF and IMFINSE have same number of columns as IMF. If
%   IMF is a timetable, IMFINSF and IMFINSE are also timetable.
%
%   [P,F,T,IMFINSF,IMFINSE] = HHT(..., 'Name1', Value1, 'Name2', Value2,
%   ...) specifies additional properties as name-value pairs. The supported
%   Name-Value pairs are:
%
%       'FrequencyLimits':      Specify frequency range (unit: Hz) for
%                               computing Hilbert Spectrum. The default
%                               range is [0, Fs/2], where Fs represents the
%                               sampling frequency.
%
%       'FrequencyResolution':  Specify frequency resolution (unit: Hz) to
%                               discretize frequency range. If
%                               'FreqResolution' is not specified, the
%                               frequency range [f_low, f_high] will be
%                               divided into 100 segments, giving a
%                               resolution of (f_high-f_low)/100.
%
%       'MinThreshold':         Set elements of P to 0 when the
%                               corresponding elements of 10*log10(P) are
%                               less than MinThreshold. It also affects the
%                               output plot. The default value is -Inf.
%
%   HHT(...) with no output arguments plots the Hilbert Spectrum in the
%   current figure.
%
%   HHT(...,FREQLOCATION) controls where MATLAB displays the frequency axis
%   on the plot. This string can be either 'xaxis' or 'yaxis'. The default
%   is 'yaxis' which displays the frequency on the y-axis.
%
%   See also EMD.
%
%   %EXAMPLE 1: Visualize Hilbert spectrum of non-stationary signals
%   fs = 1000;
%   t = (0:1/fs:4)';
%   x1 = sin(2*pi*50*t) + sin(2*pi*200*t) + 0.1*randn(length(t),1);
%   x2 = sin(2*pi*50/2*t) + sin(2*pi*200/2*t) + sin(2*pi*250*t) + 0.1*randn(length(t),1);
%   x = [x1;x2];
%   IMF = emd(x);
%   hht(IMF,fs);
%

% Copyright 2017-2018 The MathWorks, Inc.

%#codegen
narginchk(1,Inf);
nargoutchk(0,5);

[IMF,fs,T,TD,F,FRange,FResol,MinThres,Method,FreqLoc,isTT,isNF] = parseAndValidateInputs(IMF, varargin{:});
isInMATLAB = coder.target('MATLAB');

% store instantaneous frequency and energy
freqIdx = zeros(size(IMF,2), length(T));
insf = zeros(length(T),size(IMF,2));
inse = zeros(length(T),size(IMF,2));

for i = 1:size(IMF,2)
    switch Method
        % for future extension
        case 'HT'
            sig = hilbert(IMF(:,i));
            energy = abs(sig).^2;
            phaseAngle = angle(sig);
    end
    
    % compute instantaneous frequency using phase angle
    omega = gradient(unwrap(phaseAngle));
    
    % convert to Hz
    omega = fs/(2*pi)*omega;
    
    % find out index of the frequency
    omegaIdx = floor((omega-F(1))/FResol)+1;
    freqIdx(i,:) = omegaIdx(:,1)';
    
    % generate distribution
    insf(:,i) = omega;
    inse(:,i) = energy;
end

% filter out points not in the frequency range
idxKeep = (freqIdx>=1) & (freqIdx<=length(F));
timeIdx = repmat(1:length(T),size(IMF,2),1);
inseFilt = inse';%inse(:);

% store energy in sparse matrix
P = sparse(freqIdx(idxKeep),timeIdx(idxKeep),inseFilt(idxKeep),length(F),length(T));
P(10*log10(P)<MinThres) = 0;

if(isTT)
    T = TD;
end

if(nargout==0 && isInMATLAB)
    hhtPlot(insf, inse, T, FRange, MinThres, FreqLoc, isNF);
end

if(isTT && isInMATLAB)
    insf = array2timetable(insf,'RowTimes',T);
    inse = array2timetable(inse,'RowTimes',T);
end

if nargout > 0
    varargout{1} = P;
end

if nargout > 1
    varargout{2} = F;
end

if nargout > 2
    varargout{3} = T;
end

if nargout > 3
    varargout{4} = insf;
end

if nargout > 4
    varargout{5} = inse;
end

end

%--------------------------------------------------------------------------
function [IMF,fs,T,TD,F,FRange,FResol,MinThres,Method,FreqLoc,isTT,isNF] = parseAndValidateInputs(IMF, varargin)
% input type checking
validateattributes(IMF,{'single','double','timetable'},{'2d','nonempty'},'hht','IMF');
isInMATLAB = coder.target('MATLAB');
isTT = isa(IMF,'timetable');

if isvector(IMF) && (~isTT)
    IMF = IMF(:);
end

if(size(IMF,1)<2)
    error(message('shared_signalwavelet:hht:general:notEnoughRows','IMF',1));
end


% check if Fs/Frequency location exist or not
fs = 2*pi;
FreqLoc = 'yaxis';
isNF = true;    % if it is normalized frequency
initVarargin = 1;
finalVarargin = length(varargin);

if(~isempty(varargin))
    if(~ischar(varargin{1}) && ~isstring(varargin{1}))
        fs = varargin{1};
        isNF = false;
        initVarargin = 2;
    end
    
    if((ischar(varargin{end}) || isstring(varargin{end}))...
            && mod(finalVarargin-initVarargin+1,2)==1)
        FreqLoc = varargin{end};
        finalVarargin = length(varargin)-1;
    end
end
validateattributes(fs,{'numeric'},{'nonnan','finite',...
    'positive','scalar'},'hht','fs');
validatestring(FreqLoc,{'yaxis','xaxis'},'hht','freqloc');

% handle timetable
if(isTT)
    signalwavelet.internal.util.utilValidateattributesTimetable(IMF, {'regular','sorted','multichannel'}, 'hht','IMF');
    [IMF, T, TD] = signalwavelet.internal.util.utilParseTimetable(IMF);
    validateattributes(T, {'single','double'},{'nonnan','finite','real'},'hht','T');
    
    % validate input frequency coincides with timetable
    if(~isNF)
        if(abs((T(2)-T(1))-1/fs)>eps)
            error(message('shared_signalwavelet:hht:general:notMatchedFreqTimetable','IMF'));
        end
    else
        fs = 1/(T(2)-T(1));
        isNF = false;
    end
else
    TD = [];
    T = (0:(size(IMF,1)-1))'/fs;
end

% data integrity checking
validateattributes(IMF,{'single','double'},{'real','finite','nonnan','nonsparse'},'hht','IMF');

% cast to double due to sparse matrix constraints
IMF = double(IMF);
T = double(T);

% parse and validate name-value pairs
defaultFRange = [];
defaultFResol = [];
defaultMinThres = -inf;
defaultMethod = 'HT';
if(isInMATLAB)
    p = inputParser;
    addParameter(p,'FrequencyLimits',defaultFRange);
    addParameter(p,'FrequencyResolution',defaultFResol);
    addParameter(p,'MinThreshold',defaultMinThres);
    addParameter(p,'Method',defaultMethod);
    parse(p,varargin{initVarargin:finalVarargin});
    FRange = p.Results.FrequencyLimits;
    FResol = p.Results.FrequencyResolution;
    MinThres = p.Results.MinThreshold;
    Method = p.Results.Method;
else
    coder.varsize('FRange',2);
    coder.varsize('FResol');
    parms = struct( 'FrequencyLimits',           uint32(0), ...
                    'FrequencyResolution',      uint32(0), ...
                    'MinThreshold',             uint32(0), ...
                    'Method',                   uint32(0));
    pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
    FRange = eml_get_parameter_value(pstruct.FrequencyLimits,defaultFRange,varargin{initVarargin:finalVarargin});
    FResol = eml_get_parameter_value(pstruct.FrequencyResolution,defaultFResol,varargin{initVarargin:finalVarargin});
    MinThres = eml_get_parameter_value(pstruct.MinThreshold,defaultMinThres,varargin{initVarargin:finalVarargin});
    Method = eml_get_parameter_value(pstruct.Method,defaultMethod,varargin{initVarargin:finalVarargin});
end
validateattributes(MinThres,{'numeric'},{'nonnan','scalar'},'hht','MinThreshold');
validatestring(Method,{'HT','DQ'},'hht','Method');

% compute frequency range and resolution when they are not specified
if(isempty(FRange))
    FRange = [0;fs/2];
end
validateattributes(FRange,{'numeric'},...
    {'nonnan','finite','numel',2,'>=',0,'<=',fs/2},...
    'hht','FrequencyLimits');
if(FRange(1)>=FRange(2))
    error(message('shared_signalwavelet:hht:general:invalidFreqRange', 'FrequencyLimits'));
end

if(isempty(FResol))
    FResol = (FRange(2)-FRange(1))/100;
end
validateattributes(FResol,{'numeric'},{'nonnan','finite','scalar','>',0},'hht','FrequencyResolution');

% set up frequency vector
F = (FRange(1):FResol:FRange(2))';
end
