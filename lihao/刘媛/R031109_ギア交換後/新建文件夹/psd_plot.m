function psd_plot(ts,ys,col_solid)
%figure指定すること
%        [pxx,w] = pwelch(ys,[],[],[],1/ts,'power');
[pxx,w] = pwelch(ys,[],[],[],1/ts,'psd');
plot(w,pxx,col_solid), hold on
xlabel('Frequency (Hz)'), ylabel('Power spectral density')
%         if i == 3
%             legend('N1','N2','N3')
%         end
%        if i == 1
PSD_res_dis = sprintf('周波数分解能(PSD)： %5.2f [Hz] \n',1/ts/length(pxx)/2);
disp(PSD_res_dis)
%        end
% ピーク周波数
fs = 1/ts/length(pxx)*(0:length(pxx)-1); % Sampling freq.
[PSD_max2, PSD_index2] = max(pxx);
PSD_freq2 = fs(PSD_index2)/2;
PSD_dis2 = sprintf('ピーク周波数(PSD)： %5.2f [Hz], ピーク値： %5.2f \n',PSD_freq2,PSD_max2);
disp(PSD_dis2)
end