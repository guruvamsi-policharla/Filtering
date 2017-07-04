function [WTdat] = wtAtf(sig, fs, f1, f0)
% Author: Aleksandra Pidde a.pidde@lancaster.ac.uk
% based on Dmytro Iatsenko's code
%
% TO DO: 
% cone of influence
%
% calculating Morlet wavelet at frequency f1
% INPUT:
% sig: signal
% fs: sampling frequency 
% f1: frequency of interest
% f0: mother frequency of Morlet wavelet

N = length(sig);
Nq = ceil((N + 1) / 2); 
ff = [(0 : Nq - 1), -fliplr(1 : N - Nq)] * fs / N; 
ff = ff(:);

om = 2 * pi * f0; 
fwt = @(xi)(exp(-(1/2) * (om - xi).^2) - exp(-(1/2) * (om^2 + xi.^2))); % fft of Morlet wavelet

fx = fft(sig, N); 
fx(ff <= 0) = 0; 
p = 1;

freqw = ff * om / (2 * pi * f1); %frequencies for the wavelet function
fw = conj(fwt(2 * pi * freqw)); 
nidx = find(isnan(fw) | ~isfinite(fw));
if ~isempty(nidx) %to avoid NaNs due to numerics, e.g. sin(0)/0
    fw(nidx) = conj(fwt(2 * pi * freqwf(nidx) + 1e-14));
    nidx = find(isnan(fw) | ~isfinite(fw)); 
    fw(nidx) = 0;
end
conv = fx' .* fw(:); %convolution in the frequency domain
out = ((om / (2 * pi * f1))^(1 - p)) * ifft(conv, N); % calculate WT at each time
WTdat = out(end:-1:1);
end

