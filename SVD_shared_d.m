ms1 = XZ108_ms;
ms2 = XZ117_ms;
tstamp1 = XZ108_tstamp;
tstamp2 = XZ117_tstamp;
btimestamp;
[B,~] = BehavStruExtract('XZ124_XZ105_toy.txt', 1);
%%
avgtrace1 = mean(zscore(ms1.FiltTraces(:, logical(ms1.cell_label))),2);
avgtrace2 = mean(zscore(ms2.FiltTraces(:, logical(ms2.cell_label))),2);
f = figure;
f.Position = [100,100,1200,200];
a = axes;
a.NextPlot = 'add';
a.XLim = [0,round(max(tstamp1),-2)];
a.YLim = [min([avgtrace1;avgtrace2])-0.1, max([avgtrace1;avgtrace2])+0.1]
a.XLabel.String = 'Time(s)';
a.YLabel.String = 'dF';
AddBehavBlocks(a,B,btimestamp,{'human_interfere'},[0.6 0.6 0.6]);
plot(a,tstamp1, smoothdata(avgtrace1,1,'movmean',15), 'r');
plot(a,tstamp2, smoothdata(avgtrace2,1,'movmean',15), 'b');
saveas(f,'avgtrace.eps','epsc');

%% svd
close all
pcn = 1;
[~, m2tom1,td] = TStampAlign(tstamp1, tstamp2);
ms1_ca = zscore(ms1.FiltTraces(:, logical(ms1.cell_label)));
ms2_ca = zscore(ms2.FiltTraces(m2tom1, logical(ms2.cell_label)));
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
a.XLim = [0,round(max(tstamp1),-2)];
a.YLim = [min([ms1_ca_u(:,pcn);ms2_ca_v(:,pcn)])-0.1, max([ms1_ca_u(:,pcn);ms2_ca_v(:,pcn)])+0.1];
a.XLabel.String = 'Time(s)';
a.YLabel.String = 'AU';
AddBehavBlocks(a,B,btimestamp,{'human_interfere'},[0.6 0.6 0.6]);
plot(a,tstamp1, smoothdata(ms1_ca_u(:,pcn),1,'movmean',15), 'r');
plot(a,tstamp2(m2tom1), smoothdata(ms2_ca_v(:,pcn),1,'movmean',15), 'b');
corr(ms2_ca_v(:,pcn),ms1_ca_u(:,pcn))
saveas(f,['PC',num2str(pcn),'.eps'],'epsc');
explained1 = diag(ms1_ca_u'*ms1_ca_u)/sum(diag(ms1_ca_u'*ms1_ca_u));
explained2 = diag(ms2_ca_v'*ms2_ca_v)/sum(diag(ms2_ca_v'*ms2_ca_v));
f2 = figure; plot(explained1,'r'), hold on, plot(explained2,'b');
xlabel('PC#'); ylabel('Fraction of total variance explained');
saveas(f2, 'var_explained.eps','epsc');
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