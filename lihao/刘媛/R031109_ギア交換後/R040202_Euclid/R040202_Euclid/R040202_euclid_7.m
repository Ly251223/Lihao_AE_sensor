%for ii=1:1
    clear;clc;close all

    
        % 各種プロット、データ保存のON/OFF
    % 個別の軸のプロット
    show_each_axis = 0;
    % 全軸の接続プロット
    show_all_axis = 1;
    % ショットの開始同期
    sig_proc = 1;
    % PSDプロット
    show_psd = 0;
    % データ分類
    get_file = 1;9:00

    % プロットするショット数
    n_shot = 3; % max3ショット
    % 軸数
    n_axis = 4;
    % デシメーション有無
    decim = 0;
    deci = 1; % デシメーション次数
    % パス設定
    filepath = pwd;

    % データの読み込み
%     filename = '97XYUV';
%     sheet = 1;
%     xlRange = 'A1:GR2144'; % sheet = 1
%     dat = xlsread(filename,sheet,xlRange);
%     dat = dat(:,1:100); % ミラーデータの削除（100列以上は解析不要）

for n_shot = 1:3

    % データの読み込み
    cd .\dat 
    if n_shot == 1
        filename = '20190327_120441_009_00008.csv';% 8ショット目
    elseif n_shot == 2
        filename = '20190327_123539_009_00332.csv';% 332ショット目
    elseif n_shot ==3
        filename = '20190327_124520_009_00438.csv';% 438ショット目
    end
    dat = csvread(filename,4,1);
    cd ..

    %% 軸データ作成
    % サンプリング時間
    Ts = 1e-3; % 1ms
    % データ長
    length_dat = 0.4/Ts;
    %

    dat_X = dat(:,2)'; % X軸の1ショット分×行数 [x0.01MPa]
    dat_Y = dat(:,3)'; % Y軸の1ショット分×行数
    dat_U = dat(:,4)'; % U軸の1ショット分×行数
    dat_V = dat(:,5)'; % V軸の1ショット分×行数
    dat_X_cmd = dat(:,6)'; % X軸指令の1ショット分×行数
    dat_X_pos = dat(:,7)'; % D/C位置[um]
    dat_B_pos = dat(:,8)'; % スライド位置[um]
    dat_T_cmd = dat(:,9)'*0.1; % D/C上流右モータトルクコマンド[%] D/C上游右扭矩指令
    
    %% スタート同期
    if sig_proc == 1
        diff_dat_X = diff(dat_X);

        n_start = min(find(diff_dat_X>4))-3;
        n_end = length(dat_X);
        
        dat_X = dat_X(1,n_start:n_end);
        dat_Y = dat_Y(1,n_start:n_end);
        dat_U = dat_U(1,n_start:n_end);
        dat_V = dat_V(1,n_start:n_end);
        dat_X_cmd = dat_X_cmd(1,n_start:n_end);
        dat_X_pos = dat_X_pos(1,n_start:n_end);
        dat_B_pos = dat_B_pos(1,n_start:n_end);
        dat_T_cmd = dat_T_cmd(1,n_start:n_end);
    end

    % 立ち上がりの切り出し
    dat_X = dat_X(1,1:length_dat); % X軸の1ショット分×行数 [x0.01MPa]
    dat_Y = dat_Y(1,1:length_dat); % Y軸の1ショット分×行数
    dat_U = dat_U(1,1:length_dat); % U軸の1ショット分×行数
    dat_V = dat_V(1,1:length_dat); % V軸の1ショット分×行数
%    dat_X_cmd = dat(:,6); % X軸指令の1ショット分×行数
%    dat_X_pos = dat(:,7); % D/C位置[um]
%    dat_B_pos = dat(:,8); % スライド位置[um]
%    dat_T_cmd = dat(:,9)*0.1; % D/C上流右モータトルクコマンド[%]

    %% 単位調整
    dat_X = dat_X*0.01; % [x0.01MPa]⇒[MPa]
    dat_Y = dat_Y*0.01;
    dat_U = dat_U*0.01;
    dat_V = dat_V*0.01;

    %% 軸データ作成
    dat_XYUV(n_shot,:) = [dat_X(:,:),dat_Y(:,:),dat_U(:,:),dat_V(:,:)]; % 行：ショット数，列：XYUV

end    


    %% ユークリッド距離の算出
    % 全軸結合データでのショット間差
    [nc,nr]=size(dat_X);
    %tim = [0:Ts:Ts*(nr-1)];
    num = [1:1:nr];
    n_XYUV = length([dat_X(1,:),dat_Y(1,:),dat_U(1,:),dat_V(1,:)]);
    num_XYUV = [1:1:n_XYUV];
    % 軸別のショット間距離差
    % 距離
    for k=1:n_axis
        for i=1:n_shot
            dist_Euc_axis(1,i) = pdist2(dat_XYUV(1,1+400*(k-1):400+400*(k-1)),dat_XYUV(i,1+400*(k-1):400+400*(k-1)),'euclidean');
        end
        figure
        for i=1:3
            plot([0:Ts:(400-1)*Ts],dat_XYUV(i,1+400*(k-1):400+400*(k-1))),hold on
            xlabel('Time [s]'), ylabel('Amplitude [MPa]')
           % xlabel('Number of samples'), ylabel('Amplitude: X axis [MPa]')
        end
        legend('8th','332nd','438th')
        if k==1
            title('X axis')
        elseif k==2
            title('Y axis')
        elseif k==3
            title('U axis')
        elseif k==4
            title('V axis')
        end
        figure
        plot([8, 332, 438],dist_Euc_axis,'o-')
        xlabel('Number of shot'), ylabel('Euclidian distance')
        if k==1
            title('X axis')
        elseif k==2
            title('Y axis')
        elseif k==3
            title('U axis')
        elseif k==4
            title('V axis')
        end
    end

    
    
    
    % 4軸同時プロット（順次接続）
    if show_all_axis == 1
        figure
        for i=1:n_shot
            if i<=1
                plot(num_XYUV, dat_XYUV(i,:),'b'),hold on
            elseif i<=2
                plot(num_XYUV, dat_XYUV(i,:),'g')
            elseif i<=3
                plot(num_XYUV, dat_XYUV(i,:),'m')
            else
                plot(num_XYUV, dat_XYUV(i,:),'r')
            end
        end
        xlabel('Number of samples:X+Y+U+V'), ylabel('Amplitude: each shot [MPa]')
        grid on
    end
    legend('8th','332nd','438th')
    % 距離
    for i=1:n_shot
        dist_Euc(1,i) = pdist2(dat_XYUV(1,:),dat_XYUV(i,:),'euclidean');
    end
    figure
    plot([8, 332, 438],dist_Euc,'o-')
    xlabel('Number of shot'), ylabel('Euclidian distance')

    % 同じショットでの軸間差
    % ショット同士のプロット
    for k=1:3 % ショット
        figure
        for i=1:4
            plot(dat_XYUV(k,1+(400*(i-1)):400+(400*(i-1)))),hold on
            xlabel('Number of samples'), ylabel('Amplitude: each axis [MPa]')
            if k==1
                title('8th shot')
            elseif k==2
                title('332nd shot')
            elseif k==3
                title('438th shot')
            end
        end
        legend('X','Y','U','V')
    end
    % 距離
    for k=1:n_shot
        for i=1:n_axis
            dist_Euc(k,i) = pdist2(dat_XYUV(k,1:400),dat_XYUV(k,1+(400*(i-1)):400+(400*(i-1))),'euclidean');
        end
        figure
        plot(dist_Euc(k,:),'o-')
        xlabel('Number of axis'), ylabel('Euclidian distance to X-axis')
        if k==1
            title('8th shot')
        elseif k==2
            title('332nd shot')
        elseif k==3
            title('438th shot')
        end
        figure
        X = categorical({'X','Y','U','V'});
        X = reordercats(X,{'X','Y','U','V'});
        bar(X, dist_Euc(k,:))
        xlabel('Axis'), ylabel('Euclidian distance to X-axis')
        if k==1
            title('8th shot')
        elseif k==2
            title('332nd shot')
        elseif k==3
            title('438th shot')
        end
    end
    
   
    
    

    
    
%end
    
    
