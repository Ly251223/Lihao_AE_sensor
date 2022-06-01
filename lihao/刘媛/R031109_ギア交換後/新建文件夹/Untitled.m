sigLen = 4096;
fb60 = cwtfilterbank('SignalLength',sigLen);
fb10 = cwtfilterbank('SignalLength',sigLen,'TimeBandwidth',10);
% Obtain the time-domain wavelets for the filter banks.

[psi60,t] = wavelets(fb60);
[psi10,~] = wavelets(fb10);
% Use the scales function to find the mother wavelet for each filter bank.

sca60 = scales(fb60);
sca10 = scales(fb10);
[~,idx60] = min(abs(sca60-1));
[~,idx10] = min(abs(sca10-1));
m60 = psi60(idx60,:);
m10 = psi10(idx10,:);
% Since the time-bandwidth product is larger for the fb60 filter bank, 
% verify the m60 wavelet has more oscillations under its envelope than the m10 wavelet.

subplot(2,1,1)
plot(t,abs(m60))
grid on
hold on
plot(t,real(m60))
plot(t,imag(m60))
xlim([-30 30])
legend('abs(m60)','real(m60)','imag(m60)')
title('TimeBandwidth = 60')
subplot(2,1,2)
plot(t,abs(m10))
grid on
hold on
plot(t,real(m10))
plot(t,imag(m10))
xlim([-30 30])
legend('abs(m10)','real(m10)','imag(m10)')
title('TimeBandwidth = 10')