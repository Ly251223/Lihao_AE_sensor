for ii=1:1
    clear;clc;%close all

    
        % 各種プロット、データ保存のON/OFF
    % 異常値とゼロ飛びの補正
    sig_proc = 0;
    % 個別の軸のプロット
    show_each_axis = 0;
    % 全軸の接続プロット
    show_all_axis = 1;
    % PSDプロット
    show_psd = 0;
    % データ分類
    get_file = 1;
    % プロットするショット数
    n_shot = 2144; % max2144ショット
    % デシメーション有無
    decim = 0;
    deci = 1; % デシメーション次数
    % パス設定
    filepath = pwd;

    % データの読み込み
    filename = '97XYUV';
    sheet = 1;
    xlRange = 'A1:GR2144'; % sheet = 1
    dat = xlsread(filename,sheet,xlRange);
    dat = dat(:,1:100); % ミラーデータの削除（100列以上は解析不要）

    %% 異常値とゼロ飛びの補正
    if sig_proc == 1
        %% 列の削除（すべてのショット値が異常）
         dat(276,:) = []; % 毎ショット異常
         dat(517-1,:) = []; % 4ショット目で異常値

        dat(1016-2,:) = []; % ゼロ出力
         dat(1976-3,:) = []; % ゼロ出力

        %% 飛びの補正
        [nc_d,nr_d]=size(dat);
        n_mod = 0; % 補正数カウント
        for i =1:nc_d
            for j=2:nr_d
                if dat(i,j)==0 % ゼロのとき ⇒ 一つ前のデータを保持
        %        if dat(i,j)<1e-42 % 1e-42以下のとき ⇒ 一つ前のデータを保持
                    dat(i,j) = dat(i,j-1);
                    n_mod = n_mod+1;
                end
            end
        end
            Xdis1 = sprintf('データ補正数 %d [個] \n',n_mod);
            disp(Xdis1)
    end
    
    %% 小数点を移動
    dat_e38 = dat*1e38; % 小数点シフト

    %% 軸データ作成
    dat_X = dat_e38(:,1:25); % X軸の1ショット分×行数
    dat_Y = dat_e38(:,26:50); % Y軸の1ショット分×行数
    dat_U = dat_e38(:,51:75); % U軸の1ショット分×行数
    dat_V = dat_e38(:,76:100); % V軸の1ショット分×行数
    dat_XYUV = [dat_X(:,:),dat_Y(:,:),dat_U(:,:),dat_V(:,:)];

    [nc,nr]=size(dat_X);
    %tim = [0:Ts:Ts*(nr-1)];
    num = [1:1:nr];
    n_XYUV = length([dat_X(1,:),dat_Y(1,:),dat_U(1,:),dat_V(1,:)]);
    num_XYUV = [1:1:n_XYUV];

    %% 4軸同時プロット（順次接続）
    if show_all_axis == 1
        figure
        for i=1:n_shot
            if i<=500
                plot(num_XYUV, dat_XYUV(i,:),'b'),hold on
            elseif i<=1000
                plot(num_XYUV, dat_XYUV(i,:),'g')
            elseif i<=1500
                plot(num_XYUV, dat_XYUV(i,:),'m')
            else
                plot(num_XYUV, dat_XYUV(i,:),'r')
            end
        end
        xlabel('Number of samples:X+Y+U+V'), ylabel('Amplitude: each axis [*e38]')
        grid on
    end
    
    %% データ保存の有無
    % 4軸同時で異常（ゼロ出力，ゼロ飛び，1点異常）の分類
    if get_file == 1
        %dat_axis = dat_V;
        %num_axis = num;
        dat_axis = dat_XYUV;
        num_axis = num_XYUV;
        [nr_axis,nc_axis] = size(dat_axis);
        flag_diag = zeros(nr_axis,1);
        count_NG = 0; count_OK = 0;
        % pngファイル化
        figure
        %for i=1:100:length(dat_X)
        for i=1:nr_axis
            plot(num_axis, dat_axis(i,:))
            %axis([0 100 1 1.5])
            axis([0 100 0 1.5])
            for j=1:nc_axis
                %if dat_axis(i,j)==0 % ゼロのとき ⇒ 一つ前のデータを保持
             %   if dat_axis(i,j)<3.6e-6 % 以下のとき ⇒ 一つ前のデータを保持
                if (dat_axis(i,j)<3.6e-6)|(dat_axis(i,84)<0.12)|(dat_axis(i,87)<0.12) % 以下のとき ⇒ 一つ前のデータを保持
                    flag_diag(i,1) = 1;
                   % i,j
                end
            end
            % 異常値(flag_daig==1)
            if flag_diag(i,1) == 1
                count_NG = count_NG + 1;
                dat_NG(count_NG,:) = dat_axis(i,:);
            else
                count_OK = count_OK + 1;
                dat_OK(count_OK,:) = dat_axis(i,:);
            end
        end
        %xlabel('Number of samples'), ylabel('Amplitude: X axis [*e38]')
    end
    dis1 = sprintf('OKデータ数： %d [個]，NGデータ数： %d [個] \n',nr_axis-count_NG,count_NG);
    disp(dis1)   
    
    figure
    plot(dat_NG')
    xlabel('Number of samples'),ylabel('Amplitude of NG data')
    figure
    plot(dat_OK')
    xlabel('Number of samples'),ylabel('Amplitude on OK data')
    
    
    
    %% 学習データと検証データの分離
    num_teach = 1600; % 学習用データ数（1701以内）。残りの正常＋全異常が検証用
    dat_teach = dat_OK(1:num_teach,:);
    dat_valid = vertcat(dat_OK(num_teach+1:length(dat_OK),:),dat_NG);
    
    % 学習用の特徴量作成
    feature = dat_teach'; % (特徴量100点:4軸 x ショット数)

   
    figure
    plot(feature)
    xlabel('Number of selected feature'),ylabel('Values of feature')

    %% Image clustering using k-means after feature extraction with Darknet-19
    Labels(1:num_teach,1) = categorical(["1_OK"]); % ラベル作成
    numClass = 1; % クラス数=1(正常のみ)
    %C=kmeans(feature',numClass,"Start","plus");
    %[C,C2] = kmeans(feature',numClass,"Start","plus");
    [C,C2] = kmeans(feature',numClass,'Start','plus');



    %% 検証データの作成 -------------------------------------------------------------
    % 検証用の特徴量抽出
    feature_val = dat_valid'; % 残りのOKデータ＋全NGデータ
    % ラベル付け
    Labels_val(1:length(dat_OK)-num_teach,1)=categorical(["1_OK"]);
    Labels_val(length(dat_OK)-num_teach+1:length(dat_NG)+length(dat_OK)-num_teach,1)=categorical(["2_NG"]);
    
    figure
    plot(feature_val)
    xlabel('Number of selected feature:val'),ylabel('Values of selected feature:val')
   
    % 既存のクラスターを使用して検定データを分類
    % 各検定データ点から最も近い重心を求める
    % Y 内の各観測値について X 内の観測値に対するペアワイズ距離を最小のものから K 個、昇順で返す
    % D = pdist2(X,Y,Distance,'Smallest',K) 
    % C2:重心, feature_val':検証データの特徴量, val_test:距離（近い方）
    % 近い重心からの距離
    [val_test_b2,idx_test_b2] = pdist2(C2,feature_val','euclidean','Smallest',1); 
    
    
%% 距離の閾値で2分類化する    (ユークリッド距離)
    % 距離を確認する（正規化なし）
    figure
    for i=1:length(Labels_val)
        if Labels_val(i) == '1_OK'
            p1 = plot(i,val_test_b2(1,i),'ok');hold on
%         if label_count == 0;
%             plot(i,d2_mahal(i),'ob','DisplayName','none');hold on
%             label_count = 1;
%         end
        elseif Labels_val(i) == '2_NG'
            p2 = plot(i,val_test_b2(1,i),'or');hold on
        end
    end
    legend([p1 p2],{'OK','NG'})   
    xlabel('Number of data'),ylabel('Euclidian distance')
    % 正規化プロット
    figure
    for i=1:length(Labels_val)
        if Labels_val(i) == '1_OK'
            p1 = plot(i,val_test_b2(1,i)/norm(C2(1,:)),'ok');hold on
        elseif Labels_val(i) == '2_NG'
            p2 = plot(i,val_test_b2(1,i)/norm(C2(1,:)),'or');hold on
        end
    end
    legend([p1 p2],{'OK','NG'})   
    xlabel('Number of data'),ylabel('Euclidian distance')
    thresh = input('Type the threshold: ')
    % 商品有無での距離計算：
    disp('ユークリッド距離')
    no_item_max = max(val_test_b2(1,1:length(dat_OK)-num_teach)./norm(C2(1,:)));
    item_min = min(val_test_b2(1,length(dat_OK)-num_teach+1:length(val_test_b2))./norm(C2(1,:)));
    temp_x3 = ['　商品有無での距離（マイナス値は重なり有り）： ',num2str(item_min-no_item_max)];
    disp(temp_x3)

    % d_Euclid_max = maxk(val_test_b2./norm(C2(1,:)),1);% 1番目に大きい数
    d_Euclid_max = maxk(val_test_b2./norm(C2(1,:)),4); d_Euclid_max = d_Euclid_max(4); % 4番目に大きい数
    temp_x32 = ['　最大距離で正規化した距離： ',num2str( num2str( (item_min-no_item_max)/d_Euclid_max ) )];
    disp(temp_x32)
    
    %     
    kk = 1;kkk = 1;
    for j=1:length(idx_test_b2)
        if (val_test_b2(1,j)/norm(C2(1,:))) < thresh
            val_test_b2_c1(1,kk) = val_test_b2(1,j)/norm(C2(1,:)); % 学習画像からの距離で正規化
            idx_test_b2(1,j) = 1;
            kk = kk + 1;
        else
            val_test_b2_c2(1,kkk) = val_test_b2(1,j)/norm(C2(1,:)); % 学習画像からの距離で正規化
            idx_test_b2(1,j) = 2;
            kkk = kkk + 1;
        end
    end
    

    
    %% Mahalanobis distance
    d2_mahal_sq = mahal(feature_val',feature'); % 検証データ分(feature_valの行数分)の距離が算出される
    d2_mahal = sqrt(d2_mahal_sq); % ルート計算
    figure
%    label_count = 0;
    for i=1:length(Labels_val)
        if Labels_val(i) == '1_OK'
            p1 = plot(i,d2_mahal(i),'ok');hold on
%         if label_count == 0;
%             plot(i,d2_mahal(i),'ob','DisplayName','none');hold on
%             label_count = 1;
%         end
        elseif Labels_val(i) == '2_NG'
            p2 = plot(i,d2_mahal(i),'or');hold on
        end
    end
    legend([p1 p2],{'OK','NG'})   
    xlabel('Number of data'),ylabel('Mahalanobis distance')
   % legend(Labels_val(1),Labels_val(2))
   
    % 商品有無での距離計算：
    disp('マハラノビス距離')
    no_item_max_mahal = max(d2_mahal(1:length(dat_OK)-num_teach,1));
    item_min_mahal = min(d2_mahal(length(dat_OK)-num_teach+1:length(d2_mahal),1));
    temp_x4 = ['　商品有無での距離（マイナス値は重なり有り）： ',num2str(item_min_mahal-no_item_max_mahal)];
    disp(temp_x4)
    temp_x42 = ['　最大距離で正規化した距離： ',num2str( (item_min_mahal-no_item_max_mahal)/max(d2_mahal) )];
    disp(temp_x42)

    % MTSでのマハラノビス距離（正常で1前後になるようD^2/Mと変量数でわる）
    d_MTS = d2_mahal_sq/length(feature(:,1));
    figure
%    label_count = 0;
    for i=1:length(Labels_val)
        if Labels_val(i) == '1_OK'
            p1 = plot(i,d_MTS(i),'ok');hold on
%         if label_count == 0;
%             plot(i,d_MTS(i),'ob','DisplayName','none');hold on
%             label_count = 1;
%         end
        elseif Labels_val(i) == '2_NG'
            p2 = plot(i,d_MTS(i),'or');hold on
        end
    end
    legend([p1 p2],{'OK','NG'})   
    xlabel('Number of data'),ylabel('Mahalanobis distance')
   % legend(Labels_val(1),Labels_val(2))
    disp('MTS距離')
    no_item_max_MTS = max(d_MTS(1:length(dat_OK)-num_teach,1));
    item_min_MTS = min(d_MTS(length(dat_OK)-num_teach+1:length(d_MTS),1));
    temp_x5 = ['　商品有無での距離（マイナス値は重なり有り）： ',num2str(item_min_MTS-no_item_max_MTS)];
    disp(temp_x5)
   % d_MTS_max = maxk(d_MTS);% 1番目に大きい数
    d_MTS_max = maxk(d_MTS,4); d_MTS_max = d_MTS_max(4); % 3番目に大きい数
    temp_x52 = ['　最大距離で正規化した距離： ',num2str( (item_min_MTS-no_item_max_MTS)/d_MTS_max )];
    disp(temp_x52)    
    
   
%     
 end
