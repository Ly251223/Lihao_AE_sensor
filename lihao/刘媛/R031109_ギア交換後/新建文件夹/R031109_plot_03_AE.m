clear
% データの読み込み

gear_condition = 'clutch_off';
%gear_condition = 'empty';
%gear_condition = 'press';

cur_dir = pwd;
switch gear_condition
    case 'clutch_off'
        cd .\ギア交換後AEセンサデータ\01_クラッチOFF
        filename = 'P10m210611_131026.csv';
        dat = csvread(filename,6,1);
    case 'empty'
        cd .\ギア交換後AEセンサデータ\02_空打ち
        filename = 'P10m210611_131140.csv';
       % filename = 'P10m210611_131703.csv';
        dat = csvread(filename,6,1);
    case 'press'
        cd .\ギア交換後AEセンサデータ\03_プレス
        filename = 'PA1-210611110445.csv';dat = csvread(filename,6,3);
            y = dat(:,2); % binary data(13 bit/符号付14bit)
            % 最大値リミット
            for i=1:length(y)
                if y(i,1)>4999
                    y(i,1) = 4999;
                end
            end
            dat(:,2)=y(:,1);
       % filename = 'PA1-210611120452.csv';dat = csvread(filename,6,2);
        
end
cd(cur_dir);%return to original

%size(dat);

y = dat(:,2); % binary data(13 bit/符号付14bit)
Ts = 10e-3; % sampling time
tim = [0:Ts:(length(y)-1)*Ts]';


figure
plot(tim, y)
xlabel('Time [s]'),ylabel('Output [digit]')

y_V = y*5/8192; % 13bit フルスケール=5V ⇒ resolusion 0.61mV   

figure
plot(tim, y_V)
xlabel('Time [s]'),ylabel('Output [V]')

% パワースペクトル表示 %
figure
psd_plot(Ts,y_V,'b');
title('PSD of AE signal')
%axis([0 200 0 8e-3])

% スペクトログラム表示（時間周波数）
figure
subplot(2,1,1);
plot(tim, y_V)
xlabel('Time [s]'),ylabel('Output [V]')
subplot(2,1,2);
%pspectrum(y_V,1/Ts,'spectrogram','MinThreshold',-60);
pspectrum(y_V,1/Ts,'spectrogram','OverlapPercent',99, 'Leakage',1,'MinThreshold',-60,'TimeResolution', 20e-3) % 時間分解能10ms
%pspectrum(y_V,1/Ts,'spectrogram','OverlapPercent',0, 'Leakage',1)
%pspectrum(y_V,1/Ts,'spectrogram','OverlapPercent',0, 'Leakage',1,'MinThreshold',-60) % -60dB以下のノイズをカット
%pspectrum(y_V,1/Ts,'spectrogram','OverlapPercent',0, 'Leakage',1,'MinThreshold',-60,'TimeResolution', 10e-3) % 時間分解能10ms
colorbar 'off'
figure
%pspectrum(y_V,1/Ts,'spectrogram','OverlapPercent',99, 'Leakage',1,'MinThreshold',-60,'TimeResolution', 10e-3) % 時間分解能10ms
pspectrum(y_V,1/Ts,'spectrogram','OverlapPercent',99, 'Leakage',1,'MinThreshold',-60,'TimeResolution', 20e-3) % 時間分解能10ms

figure
%pspectrum(y_V,1/Ts,'spectrogram','MinThreshold',-60,'Leakage',1,'TimeResolution', 10e-3); % 時間分解能10ms
pspectrum(y_V,1/Ts,'spectrogram','MinThreshold',-60,'Leakage',1); % 

%% スカログラムの表示：信号アナライザで解析
%signalAnalyzer(y_V,'TimeValues',tim)

%% 個別のCWTプロット
figure
% subplot(2,1,1);
% plot(tim, y_V)
% xlabel('Time [s]'),ylabel('Output [V]')
%subplot(2,1,2);
cwt(y_V, 1/Ts)

%% AE解析
% オリジナル波形
figure
subplot(211)
plot(tim,y_V)
xlabel('Time [s]'), ylabel('AE signal [V]')
title(['AE原信号'])
% 指定時間内での最大振幅
tim_length = 0.1; % 窓間隔 [s]
window_length = tim_length/Ts; % 窓の長さ
amp_window1 = zeros(length(y_V),1); % 初期値
for i=window_length+1:length(y_V)
    amp_window1(i,1) = max(y_V(i-window_length:i,1));
    %no_Hits(i,1) - no_Hits(i-window_length,1);
end
subplot(212)
plot(tim,amp_window1)
xlabel('Time [s]'),ylabel('Max. amp. per time')
title(['過去',num2str(tim_length),'[s]間での最大振幅'])
% 指定時間内0.5sでの最大振幅
tim_length = 0.5; % 窓間隔 [s]
window_length = tim_length/Ts; % 窓の長さ
amp_window1 = zeros(length(y_V),1); % 初期値
for i=window_length+1:length(y_V)
    amp_window1(i,1) = max(y_V(i-window_length:i,1));
    %no_Hits(i,1) - no_Hits(i-window_length,1);
end
figure
subplot(211)
plot(tim,amp_window1)
xlabel('Time [s]'),ylabel('Max. amp. per time')
title(['過去',num2str(tim_length),'[s]間での最大振幅'])

% 指定時間内1sでの最大振幅
tim_length = 1; % 窓間隔 [s]
window_length = tim_length/Ts; % 窓の長さ
amp_window1 = zeros(length(y_V),1); % 初期値
for i=window_length+1:length(y_V)
    amp_window1(i,1) = max(y_V(i-window_length:i,1));
    %no_Hits(i,1) - no_Hits(i-window_length,1);
end
subplot(212)
plot(tim,amp_window1)
xlabel('Time [s]'),ylabel('Max. amp. per time')
title(['過去',num2str(tim_length),'[s]間での最大振幅'])



% 発生数（閾値以上のカウントアップ）
Thresh_amp = 0.3;
Thresh_count = 0;
for i=1:length(y_V)
    if y_V(i,1) > Thresh_amp
        Thresh_count = Thresh_count +1;
    end
    no_Hits(i,1) = Thresh_count;
end
figure
subplot(211)
plot(tim,y_V)
xlabel('Time [s]'), ylabel('AE signal [V]')
title(['AE原信号'])
subplot(212)
plot(tim,no_Hits)
xlabel('Time [s]'),ylabel('Hits')
title(['発生数（閾値:',num2str(Thresh_amp),'）'])

% 過去0.1sあたりの発生数
tim_length = 0.1; % 窓間隔 [s]
window_length = tim_length/Ts; % 窓の長さ Ts = sampling time
no_Hits_window1 = zeros(length(y_V),1); % 初期値
for i=window_length+1:length(y_V)
    no_Hits_window1(i,1) = no_Hits(i,1) - no_Hits(i-window_length,1);
end
figure
subplot(211)
plot(tim,no_Hits_window1)
xlabel('Time [s]'),ylabel('Hits')
title(['過去',num2str(tim_length),'[s]間での発生数（閾値:',num2str(Thresh_amp),'）'])

% 過去0.5sあたりの発生数
tim_length = 0.5; % 窓間隔 [s]
window_length = tim_length/Ts; % 窓の長さ
no_Hits_window = zeros(length(y_V),1); % 初期値
for i=window_length+1:length(y_V)
    no_Hits_window(i,1) = no_Hits(i,1) - no_Hits(i-window_length,1);
end
subplot(212)
plot(tim,no_Hits_window)
xlabel('Time [s]'),ylabel('Hits')
title(['過去',num2str(tim_length),'[s]間での発生数（閾値:',num2str(Thresh_amp),'）'])

% 過去1sあたりの発生数
tim_length = 1; % 窓間隔 [s]
window_length = tim_length/Ts; % 窓の長さ
no_Hits_window1 = zeros(length(y_V),1); % 初期値
for i=window_length+1:length(y_V)
    no_Hits_window1(i,1) = no_Hits(i,1) - no_Hits(i-window_length,1);
end
figure
subplot(211)
plot(tim,no_Hits_window1)
xlabel('Time [s]'),ylabel('Hits')
title(['過去',num2str(tim_length),'[s]間での発生数（閾値:',num2str(Thresh_amp),'）'])

% 過去5sあたりの発生数
tim_length = 5; % 窓間隔 [s]
window_length = tim_length/Ts; % 窓の長さ
no_Hits_window = zeros(length(y_V),1); % 初期値
for i=window_length+1:length(y_V)
    no_Hits_window(i,1) = no_Hits(i,1) - no_Hits(i-window_length,1);
end
subplot(212)
plot(tim,no_Hits_window)
xlabel('Time [s]'),ylabel('Hits')
title(['過去',num2str(tim_length),'[s]間での発生数（閾値:',num2str(Thresh_amp),'）'])


% エネルギー（AI波形の面積）
area = trapz(tim,abs(y_V)) % 面積合計
cum_area = cumtrapz(tim,abs(y_V)); % 積算面積
figure
subplot(211)
plot(tim,y_V)
xlabel('Time [s]'), ylabel('AE signal [V]')
title(['AE原信号'])
subplot(212)
plot(tim,cum_area)
xlabel('Time [s]'),ylabel('AE Energy (=Area)')
title('AEエネルギー（波形の面積）')
% 過去0.1sあたりの面積
tim_length = 0.1; % 窓間隔 [s]
window_length = tim_length/Ts; % 窓の長さ
area_window = zeros(length(y_V),1); % 初期値
for i=window_length+1:length(y_V)
    area_window(i,1) = trapz(tim(i-window_length:i,1),abs(y_V(i-window_length:i,1)));
end
figure
subplot(211)
plot(tim,area_window)
xlabel('Time [s]'),ylabel('AE Energy per time')
title(['過去',num2str(tim_length),'[s]間でのAEエネルギー'])
% 過去0.5sあたりの面積
tim_length = 0.5; % 窓間隔 [s]
window_length = tim_length/Ts; % 窓の長さ
area_window = zeros(length(y_V),1); % 初期値
for i=window_length+1:length(y_V)
    area_window(i,1) = trapz(tim(i-window_length:i,1),abs(y_V(i-window_length:i,1)));
end
subplot(212)
plot(tim,area_window)
xlabel('Time [s]'),ylabel('AE Energy per time')
title(['過去',num2str(tim_length),'[s]間でのAEエネルギー'])



% 周波数
figure
subplot(211)
plot(tim,y_V)
xlabel('Time [s]'), ylabel('AE signal [V]')
title(['AE原信号'])
subplot(212)
psd_plot(Ts,y_V,'b');
title('パワースペクトル')










