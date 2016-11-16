function [ B, A ] = lpApply( Input, Order, Freq, color )
%lpApply This function plots the given datapoints after being run through
%a LPfilter.
if nargin < 4|| isempty(color)
   color = 'b-';
end

[NFFT, ~] =size(Input);
f = 1/2*linspace(0,1,NFFT/2+1);
[B,A] = butter(Order,Freq,'low');
FFTsignal = fft(Input, NFFT);
FilterFreqResponse = freqz(B,A,NFFT,'whole');
FFTfilteredSignal = FFTsignal .* FilterFreqResponse;
filteredSignal = ifft(FFTfilteredSignal, NFFT);
%filteredSignal = ifft(FFTsignal, NFFT);
x = linspace(0,NFFT,NFFT);
plot(x,filteredSignal,color);

end

