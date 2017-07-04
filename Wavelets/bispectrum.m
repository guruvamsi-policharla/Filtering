function [bxxx, bppp, bpxx, bxpp] = bispectrum(sig1, sig2, wt_1, wt_2, freqarr, eps, ord, fs, f0)
%BISPECTRUM Summary of this function goes here

bxxx = NaN*zeros(length(freqarr));
bppp = NaN*zeros(length(freqarr));
bpxx = NaN*zeros(length(freqarr));
bxpp = NaN*zeros(length(freqarr));
h = waitbar(0,'Please wait...');
tic
count = 0;
for i = 1:length(freqarr)
    for j = 1:length(freqarr)        

        f3 = freqarr(i)+freqarr(j);        
        
        if(f3 <= freqarr(end)) && abs((log2(freqarr(j)) - log2(freqarr(i))) <= ord) && abs((log2(freqarr(i)) - log2(freqarr(j))) <= ord)
            count = count + 1;
            k = find(freqarr >= f3, 1);
            if abs(freqarr(k) - f3) > eps                
                WT1f3 = wtAtf(sig1, fs, f3, f0).';
                WT2f3 = wtAtf(sig2, fs, f3, f0).';             
                bxxx(i, j) = abs(nanmean(wt_1(i, :).*wt_1(j, :).*conj(WT1f3)));
                bppp(i, j) = abs(nanmean(wt_2(i, :).*wt_2(j, :).*conj(WT2f3)));
                bxpp(i, j) = abs(nanmean(wt_1(i, :).*wt_2(j, :).*conj(WT2f3)));
                bpxx(i, j) = abs(nanmean(wt_2(i, :).*wt_1(j, :).*conj(WT1f3)));
            else
                bxxx(i, j) = abs(nanmean(wt_1(i, :).*wt_1(j, :).*conj(wt_1(k, :))));    
                bppp(i, j) = abs(nanmean(wt_2(i, :).*wt_2(j, :).*conj(wt_2(k, :))));
                bxpp(i, j) = abs(nanmean(wt_1(i, :).*wt_2(j, :).*conj(wt_2(k, :))));
                bpxx(i, j) = abs(nanmean(wt_2(i, :).*wt_1(j, :).*conj(wt_1(k, :))));              
            end
        end       
    end
    waitbar(i / length(freqarr));
end

close(h); 
toc
end

