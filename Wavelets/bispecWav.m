function Bisp = bispecWav(sig2, WT1, WT2, freq, ord, eps, f0, fs)
% Author: Aleksandra Pidde a.pidde@lancaster.ac.uk
%
% calculating wavelet bispectrum
% sig2: second signal
% WT1, WT2: wavelet transforms,
% freq: frequency axis
% eps: approximation for the third scale of the frequncy, the smaller the more
% accurate

%display = true;
display = false;
nfreq = length(freq);
Bisp = NaN * zeros(nfreq, nfreq);
auto = false;
% if compareMatrix(WT1, WT2)
%     auto = true;
% end
h = waitbar(0,'Please wait...');
tic

for j = 1 : nfreq
    kstart = 1;
    if auto
        kstart = j;
    end
    for k = kstart : nfreq
        f3 = freq(j) + freq(k);
        if (f3 <= freq(end)) && abs((log2(freq(k)) - log2(freq(j))) <= ord) && abs((log2(freq(j)) - log2(freq(k))) <= ord)
            %count = count +1;           
            idx3 = find(freq >= f3, 1);
            if abs(freq(idx3) - f3) > eps
                WTdat = wtAtf(sig2, fs, f3, f0);
                WTdat = WTdat(:)'; % make sure it is horizontal vector
                if display; disp(['calculating for a scale real ' num2str(f3) ' vs ' num2str(freq(idx3))]); end
            else
                WTdat = WT2(idx3, :);
            end
            xx = WT1(j, :) .* WT2(k, :) .* conj(WTdat);
            %xx = TFR1(j, :) .* TFR2(k, :) * transpose(conj(TFR2(idx3, :)));
            ss = nanmean(xx);
            Bisp(j, k) = abs(ss);
        end
    end
    waitbar(j / nfreq);
end
display(count1);
display(count2);
close(h); 
toc
end
