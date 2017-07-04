%%% Autor: Aleksandra Pidde aleksandra.pidde@gmail.com,
%%% a.pidde@lancaster.ac.uk

function [biamp, biphase] = biphaseWav(X, TFR1, TFR2, freq, f1, f2, fs, f0)
% Calculation of a biamplitude and biphase
% X: time series
% TFR1, TFR2: time-frequency representations 
% freq: frequency scale
% f1, f2: frequencies of interest
% fs: sampling frequency
% f0: central frequency

%fprintf('\nBisp Acommplished:  0 %');
%cont = 10;
tic
idx1 = find(freq <= f1);
idx1 = idx1(end);
idx2 = find(freq <= f2);
idx2 = idx2(end);

f3 = freq(idx1) + freq(idx2);
WTdat = wtAtf(X, fs, f3, f0);
WTdat = WTdat(:)';

xx = TFR1(idx1, :) .* TFR2(idx2, :) .* conj(WTdat);
biamp = abs(xx);
%biphase = unwrap(angle(TFR1(idx1, :))) + unwrap(angle(TFR2(idx2, :))) - unwrap(angle(conj(TFR2(idx3, :))));
biphase = unwrap(angle(xx));
end


