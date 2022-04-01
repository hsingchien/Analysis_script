this_pair = DLX{1};
ms1 = this_pair{1}.MS{2};
ms2 = this_pair{2}.MS{2};
m2tom1 = this_pair{1}.TimeStamp.mapTs{2}.M2toM1;
tstamp1 = this_pair{1}.TimeStamp.Ts{2}.Ms;
tstamp2 = this_pair{2}.TimeStamp.Ts{2}.Ms;
btimestamp = this_pair{2}.TimeStamp.Ts{2}.Bv;
B1 = this_pair{1}.Behavior{2};
B2 = this_pair{2}.Behavior{2};
m1 = this_pair{1}.MS{2};
m2 = this_pair{2}.MS{2};
%%
avgtrace1 = mean(zscore(ms1.FiltTraces(:, logical(ms1.goodCellVec))),2);
avgtrace2 = mean(zscore(ms2.FiltTraces(:, logical(ms2.goodCellVec))),2);
f = figure;
f.Position = [100,100,1200,200];
a = axes;
a.NextPlot = 'add';
a1.XLim = [0,round(max(tstamp1),-2)/1000];
a1.YLim = [min([avgtrace1;avgtrace2])-0.1, max([avgtrace1;avgtrace2])+0.1]
a1.XLabel.String = 'Time(s)';
a1.YLabel.String = 'dF';
AddBehavBlocks(a1,B1,btimestamp/1000,social_bv,assignColors(social_bv));
plot(a1,tstamp1/1000, smoothdata(avgtrace1,1,'movmean',15), 'r-');
plot(a1,tstamp2/1000, smoothdata(avgtrace2,1,'movmean',15), 'b-');
% saveas(f,'avgtrace.eps','epsc');
%% PCA
coeff1 = pca(zscore(ms1.FiltTraces(:, logical(ms1.goodCellVec))));
coeff2 = pca(zscore(ms2.FiltTraces(:, logical(ms2.goodCellVec))));
ms1_proj = zscore(ms1.FiltTraces(:, logical(ms1.goodCellVec))) * coeff1;
ms2_proj = zscore(ms2.FiltTraces(:, logical(ms2.goodCellVec))) * coeff2;
figure, a2 = axes; a2.NextPlot = 'add';
plot(a2, tstamp1/1000, ms1_proj(:,1),'r-'); plot(a2, tstamp2/1000, ms2_proj(:,1), 'b-');
AddBehavBlocks(a2,B1,btimestamp/1000,social_bv,assignColors(social_bv));
a2.XLim = [0, tstamp1(end)/1000];
corr(ms1_proj(:,1),ms2_proj(m2tom1,1));
explained1 = diag(cov(ms1_proj))/sum(diag(cov(ms1_proj)));
explained2 = diag(cov(ms2_proj))/sum(diag(cov(ms2_proj)));

%% svd
% close all
social_bv = {'attack','chasing','escape','defend','flinch','tussling','threaten',...
    'general-sniffing','sniff_face','sniff_genital',...
    'approach','follow','interaction'}; 
pcn = 1;

ms1_ca = zscore(ms1.FiltTraces(:, logical(ms1.goodCellVec)));
ms2_ca = zscore(ms2.FiltTraces(m2tom1, logical(ms2.goodCellVec)));
cov_mat = ms1_ca'*ms2_ca;
[U,S,V] = svd(cov_mat);
ms1_ca_u = ms1_ca * U;
ms2_ca_v = ms2_ca * V;
if abs(min(ms1_ca_u(:,pcn))) > abs(max(ms1_ca_u(:,pcn)))
    ms1_ca_u(:,pcn) = -1*ms1_ca_u(:,pcn);
end
if abs(min(ms2_ca_v(:,pcn))) > abs(max(ms2_ca_v(:,pcn)))
    ms2_ca_v(:,pcn) = -1*ms2_ca_v(:,pcn);
end
f = figure;
f.Position = [100,100,1200,200];
a = axes;
a.NextPlot = 'add';
a3.XLim = [0,round(max(tstamp1),-2)/1000];
a3.YLim = [min([ms1_ca_u(:,pcn);ms2_ca_v(:,pcn)])-0.1, max([ms1_ca_u(:,pcn);ms2_ca_v(:,pcn)])+0.1];
a3.XLabel.String = 'Time(s)';
a3.YLabel.String = 'AU';
AddBehavBlocks(a3,B1,btimestamp/1000,social_bv,assignColors(social_bv));
plot(a3,tstamp1/1000, smoothdata(ms1_ca_u(:,pcn),1,'movmean',15), 'r');
plot(a3,tstamp2(m2tom1)/1000, smoothdata(ms2_ca_v(:,pcn),1,'movmean',15), 'b');
corr(ms2_ca_v(:,pcn),ms1_ca_u(:,pcn))
% saveas(f,['PC',num2str(pcn),'.eps'],'epsc');
explained1 = diag(ms1_ca_u'*ms1_ca_u)/sum(diag(ms1_ca_u'*ms1_ca_u));
explained2 = diag(ms2_ca_v'*ms2_ca_v)/sum(diag(ms2_ca_v'*ms2_ca_v));
f2 = figure; plot(explained1,'r'), hold on, plot(explained2,'b');
xlabel('PLSC#'); ylabel('Fraction of total variance explained');
% saveas(f2, 'var_explained.eps','epsc');
%% permute 
ms1_pc1s = [];
ms2_pc1s = [];
for k = 1:100
    ms1_ca = zscore(ms1.FiltTraces(:, ms1.cell_label));
    ms1_ca = ms1_ca(randperm(size(ms1_ca,1)),:);
    ms2_ca = zscore(ms2.FiltTraces(m2tom1, ms2.cell_label));
%     ms2_ca = ms2_ca(randperm(size(ms2_ca,1)),:);
    cov_mat = ms1_ca'*ms2_ca;
    [U,S,V] = svd(cov_mat);
    ms1_ca_u = ms1_ca * U;
    ms2_ca_v = ms2_ca * V;
    ms1_pc1s = [ms1_pc1s, ms1_ca_u(:,1)];
    ms2_pc1s = [ms2_pc1s, ms2_ca_v(:,1)];
end
% CI95 = tinv([0.025,0.975], size(ms1_pc1s,2));
% ms1_pc1s_sem = std(ms1_pc1s,0,2)/sqrt(size(ms1_pc1s,2));
% ms2_pc1s_sem = std(ms2_pc1s,0,2)/sqrt(size(ms2_pc1s,2));
% ms1_CI95 = [ms1_pc1s_sem*CI95(1),ms1_pc1s_sem*CI95(2)];
% ms2_CI95 = [ms2_pc1s_sem*CI95(1),ms2_pc1s_sem*CI95(2)];
f = figure; a = axes; a.NextPlot = 'add';
f.Position = [100,100,1200,200];
a.XLim = [0,round(max(tstamp1),-2)];
a.YLim = [min([mean(ms1_pc1s,2);mean(ms2_pc1s,2)])-0.1, max([mean(ms1_pc1s,2);mean(ms1_pc1s,2)])+0.1];
plot(a,tstamp1, mean(ms1_pc1s,2),'r'); plot(a,tstamp2(m2tom1),mean(ms2_pc1s,2),'b');
AddBehavBlocks(a,B,btimestamp,{'human_interfere'},[0.6 0.6 0.6]);
% % patch([tstamp1;fliplr(tstamp1)],[ms1_CI95(:,1); fliplr(ms1_CI95(:,2))],'r','EdgeColor','none','FaceAlpha',0.25)
%%
%% sparseness
% percentage of off cells

cur_cell = DLX;
cur_type = 'DLX';
for i = 1:length(cur_cell)
    for j = 1:2
        cur_animal = cur_cell{i}{j};
        for m = 1:length(cur_animal.MS)
            cur_S = transpose(cur_animal.MS{m}.S);
            cur_S = cur_S(:,logical(cur_animal.MS{m}.goodCellVec));
            cur_S_threshold = mean(cur_S,1) + 2*std(cur_S);
            cur_S_filtered = double(cur_S >= repmat(cur_S_threshold,[size(cur_S,1),1]));
            twindow = 30;
            sparsenesses = [];
            pcts = [];
            tws = [];
            for tw = 1:round(twindow/2):(size(cur_S_filtered,1)-twindow+1)
                cur_tw = cur_S_filtered(tw:(tw+twindow-1),:);
                pct = sum((sum(cur_tw,1)>0))/size(cur_S_filtered,2);
                fr = sum(cur_tw,1); % 1 x n_neuron
                sparseness = (sum(fr)/size(cur_S_filtered,2))^2/sum(fr.^2/size(cur_S_filtered,2));
                tws = [tws,tw];
                pcts = [pcts, pct];
                sparsenesses = [sparsenesses, sparseness];
            end
            f = figure;
            a1 = subplot(2,1,1);
            a2 = subplot(2,1,2);
            plot(a1, tws, pcts); title(a1,'PCT');
            plot(a2, tws, sparsenesses); title(a2,'Sparseness');
            savefig(f,[cur_type,'\P',num2str(i),'A',num2str(j),'M',num2str(m)]);
            close(f);
        end
    end
end
%% get significant pairs, threshold by pcc & var 97.5%
cur_cell = CMK;
cur_type = 'CMK';
threshold = 97.5;
pair_count = 0;
ms1_passers = 0;
ms2_passers = 0;
corr_passers = 0;
parfor i = 1:length(cur_cell)
    cur_a1 = cur_cell{i}{1};
    cur_a2 = cur_cell{i}{2};
    for m = 1:length(cur_a1.MS)
        if size(cur_a1.MS{m}.FiltTraces,1) < 16000
            continue;
        end
        pair_count = pair_count+1;
        m2tom1 = cur_a1.TimeStamp.mapTs{m}.M2toM1;
        tstamp1 = cur_a1.TimeStamp.Ts{m}.Ms;
        tstamp2 = cur_a2.TimeStamp.Ts{m}.Ms;
        btimestamp = cur_a2.TimeStamp.Ts{m}.Bv;
        ms1 = cur_a1.MS{m};
        ms2 = cur_a2.MS{m};
        pcn = 1;

        ms1_ca = zscore(ms1.FiltTraces(:, logical(ms1.goodCellVec)));
        ms2_ca = zscore(ms2.FiltTraces(m2tom1, logical(ms2.goodCellVec)));
        if size(ms1_ca,1) ~= size(ms2_ca,1)
            ms1_ca = ms1_ca(1:min(size(ms1_ca,1), size(ms2_ca,1)),:);
            ms2_ca = ms2_ca(1:min(size(ms1_ca,1), size(ms2_ca,1)),:);
        end
        cov_mat = ms1_ca'*ms2_ca;
        [U,S,V] = svd(cov_mat);
        ms1_ca_u = ms1_ca * U;
        ms2_ca_v = ms2_ca * V;
        plsccor = corr(ms2_ca_v(:,pcn),ms1_ca_u(:,pcn));
        explained1 = diag(ms1_ca_u'*ms1_ca_u)/sum(diag(ms1_ca_u'*ms1_ca_u));
        explained2 = diag(ms2_ca_v'*ms2_ca_v)/sum(diag(ms2_ca_v'*ms2_ca_v));
        % permute
        ms1_pc1s = [];
        ms2_pc1s = [];
        for k = 1:1000
            % circshift ms1_ca
            shiftamount = randi(size(ms2_ca,1), 1,size(ms2_ca,2));
            ms2_ca_per = ms2_ca;
%               shiftamount = randi(size(ms1_ca,1));
            for cc = 1:size(ms2_ca,2)
                ms2_ca_per(:,cc) = circshift(ms2_ca(:,cc), shiftamount(cc));
            end
        %     ms1_ca = ms1_ca(randperm(size(ms1_ca,1)),:);
%             ms1_ca_per = circshift(ms1_ca,shiftamount,1);
        %     ms2_ca = ms2_ca(randperm(size(ms2_ca,1)),:);
            cov_mat = ms1_ca'*ms2_ca_per;
            [UU,SS,VV] = svd(cov_mat);
            ms1_ca_uu = ms1_ca * UU;
            ms2_ca_vv = ms2_ca_per * VV;
            ms1_pc1s = [ms1_pc1s, ms1_ca_uu(:,1)];
            ms2_pc1s = [ms2_pc1s, ms2_ca_vv(:,1)];
        end

        ms_corrs = diag(corr(ms1_pc1s, ms2_pc1s));
        corr_thresh = prctile(ms_corrs, threshold);
        ms1_vars = diag(cov(ms1_pc1s));
        ms1_thresh = prctile(ms1_vars, threshold);
        ms2_vars = diag(cov(ms2_pc1s));
        ms2_thresh = prctile(ms2_vars, threshold);
        if cov(ms1_ca_u(:,pcn)) >= ms1_thresh
            fprintf('Pair %d animal 1 session %d pass var\n', i, m);
            ms1_passers = ms1_passers+1;
        end
        if cov(ms2_ca_v(:,pcn)) >= ms2_thresh
            fprintf('Pair %d animal 2 session %d pass var\n', i, m);
            ms2_passers = ms2_passers+1;
        end
        if plsccor >= corr_thresh
            fprintf('Pair %d session %d pass COR\n', i, m);
            corr_passers = corr_passers + 1;
        end        
    end
end
        


      