function [S_welch] = PSD_Welch2(line,overlap,M)
%Output:Welch PSD estimate
%Input:
%1) Line:    Data in which PSD estimate is carried out
%2) Overlap: Overlap in percentage
%3) M:       Window size
D = M/100*overlap;                     %# overlapping samples
numSamples = length(line);             %Length of data
win=0.54-0.46*cos(2*pi*((1:M)/(M-1))); %M-sample Hamming
%win=rectwin(M).';
%win=hann(M).';
U=sum(win.^2)/M;                       %Window Normalization factor

%Solving for maximum index,
syms index
idx=double(floor(solve(M+index*D == numSamples,index)));

        for i=0:idx
            %Create D-sample overlapping segments, and window each segment
             seg_temp     = line(1+D*i:M+i*D).*win;
            %Compute the modified periodogram
            FFT = fft(seg_temp,M);
            PSD_mod(i+1,:) = (1/(M*U)).*(FFT.*conj(FFT)); %original
      
        end
    
%Average periodograms along first dimension to get Welch estimate
S_welch=mean(PSD_mod,1);

end
