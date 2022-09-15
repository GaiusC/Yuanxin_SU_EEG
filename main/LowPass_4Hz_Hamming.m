function Hd = LowPass_4Hz_Hamming
%LOWPASS_4HZ_HAMMING Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.11 and Signal Processing Toolbox 8.7.
% Generated on: 08-Aug-2022 02:46:47

% FIR Window Lowpass filter designed using the FIR1 function.

% All frequency values are in Hz.
Fs = 200;  % Sampling Frequency

N    = 257;      % Order
Fc   = 4;        % Cutoff Frequency
flag = 'scale';  % Sampling Flag

% Create the window vector for the design algorithm.
win = hamming(N+1);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Fc/(Fs/2), 'low', win, flag);
Hd = dfilt.dffir(b);

% [EOF]
