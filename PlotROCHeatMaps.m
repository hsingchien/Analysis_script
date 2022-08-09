% get list of mat files
flist = {dir('*.mat').name};
targetFolder = 'Tuning\';
pthreshold = 97.5;
all_active = struct(); all_active.self = {}; all_active.partner = {}; all_active.selfdrvd = {}; all_active.partnerdrvd = {};


GenTracker = {}; ATracker = {}; PairTracker = [];
%% load and plot
for i = 1:numel(flist)
    load(flist{i});
    M1 = rocStruct{1}; M2 = rocStruct{2};
    if isempty(M1.ROC)
        continue
    end
    
    if isfield(M1,'GenType')
        m1G = M1.GenType; m2G = M2.GenType;
    else
        m1G = 'C57'; m2G = 'C57';
    end
    sessionN = 'exp'; fnames = fields(M1.ROC);
    if isempty(find(contains(fnames,sessionN)))
        continue;
    end
    
    session = fnames{find(contains(fnames,sessionN))};


    if strcmp(sessionN,'exp'); expBvs = M1.ROC.(session).other.self.Behavior.EventNames; expBvs = strrep(expBvs,'_','-');end
    if strcmp(sessionN,'toy'); toyBvs = M1.ROC.(session).other.self.Behavior.EventNames; toyBvs = strrep(toyBvs,'_','-');end
    self1 = getSigCells(M1.ROC.(session).other.self.Encoding,pthreshold);
%     heatMap(self1,'xTickLabels',M1.ROC.(session).other.self.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
%     set(gcf,'Position',[0,0,800,700]); saveas(gcf,[targetFolder,sessionN,'\','self\',M1.ExperimentID,'-',M1.AnimalID,'-',m1G,'.png'],'png');
    
    self2 = getSigCells(M2.ROC.(session).other.self.Encoding,pthreshold);
%     heatMap(self2,'xTickLabels',M2.ROC.(session).other.self.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
%     set(gcf,'Position',[0,0,800,700]); saveas(gcf,[targetFolder,sessionN,'\','self\',M2.ExperimentID,'-',M2.AnimalID,'-',m2G,'.png'],'png');

    partner1 = getSigCells(M1.ROC.(session).other.partner.Encoding,pthreshold);
%     heatMap(partner1,'xTickLabels',M1.ROC.(session).other.partner.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
%     set(gcf,'Position',[0,0,800,700]); saveas(gcf,[targetFolder,sessionN,'\','partner\',M1.ExperimentID,'-',M1.AnimalID,'-',m1G,'.png'],'png');

    partner2 = getSigCells(M2.ROC.(session).other.partner.Encoding,pthreshold);
%     heatMap(partner2,'xTickLabels',M2.ROC.(session).other.partner.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
%     set(gcf,'Position',[0,0,800,700]); saveas(gcf,[targetFolder,sessionN,'\','partner\',M2.ExperimentID,'-',M2.AnimalID,'-',m2G,'.png'],'png');

    selfdrvd1 = getSigCells(M1.ROC.(session).other.derived_self.Encoding,pthreshold);
%     heatMap(selfdrvd1,'xTickLabels',M1.ROC.(session).other.derived_self.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
%     set(gcf,'Position',[0,0,600,700]); saveas(gcf,[targetFolder,sessionN,'\','self_derived\',M1.ExperimentID,'-',M1.AnimalID,'-',m1G,'.png'],'png');

    selfdrvd2 = getSigCells(M2.ROC.(session).other.derived_self.Encoding,pthreshold);
%     heatMap(selfdrvd2,'xTickLabels',M2.ROC.(session).other.derived_self.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
%     set(gcf,'Position',[0,0,600,700]); saveas(gcf,[targetFolder,sessionN,'\','self_derived\',M2.ExperimentID,'-',M2.AnimalID,'-',m2G,'.png'],'png');

    partnerdrvd1 = getSigCells(M1.ROC.(session).other.derived_partner.Encoding,pthreshold);
%     heatMap(partnerdrvd1,'xTickLabels',M1.ROC.(session).other.derived_partner.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
%     set(gcf,'Position',[0,0,600,700]); saveas(gcf,[targetFolder,sessionN,'\','partner_derived\',M1.ExperimentID,'-',M1.AnimalID,'-',m1G,'.png'],'png');
    
    partnerdrvd2 = getSigCells(M2.ROC.(session).other.derived_partner.Encoding,pthreshold);
%     heatMap(partnerdrvd2,'xTickLabels',M2.ROC.(session).other.derived_partner.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
%     set(gcf,'Position',[0,0,600,700]); saveas(gcf,[targetFolder,sessionN,'\','partner_derived\',M2.ExperimentID,'-',M2.AnimalID,'-',m2G,'.png'],'png');
    
    close all

    all_active.self = [all_active.self, self1, self2];
    all_active.partner = [all_active.partner, partner1, partner2];
    all_active.selfdrvd = [all_active.selfdrvd, selfdrvd1, selfdrvd2];
    all_active.partnerdrvd = [all_active.partnerdrvd, partnerdrvd1, partnerdrvd2];
    GenTracker = [GenTracker, m1G, m2G];
    ExID = strsplit(M1.ExperimentID,'_');
    ATracker = [ATracker, [ExID{1},'-',M1.AnimalID], [ExID{1},'-',M2.AnimalID]];
    
    PairTracker = [PairTracker, str2num(ExID{1}(regexp(ExID{1}, '\d*'):end))];
    

end

%% boxplot fraction of suppressed/active cells by genotypes
types = {'self','partner','selfdrvd','partnerdrvd'};
gens = {'KO','HET','WT'};
KOfrac = cell(1,4); HETfrac = cell(1,4); WTfrac = cell(1,4); C57frac = cell(1,4);
catBvs = {'sniff','aggressive','defensive','exploratory','initiative','investigateOBJ','other'};
for i = 1:numel(types)
   
    cur_act = all_active.(types{i});
    
    for j = 1:numel(cur_act)
        active_frac = sum(cur_act{j}==1,1)/size(cur_act{j},1);
        active_frac(isnan(cur_act{j}(1,:))) = nan; % keep nans
        suppress_frac = sum(cur_act{j}==-1,1)/size(cur_act{j},1);
        suppress_frac(isnan(cur_act{j}(1,:))) = nan;
%         if length(active_frac)==7
%             active_frac = active_frac([1:5,7]);
%             suppress_frac = suppress_frac([1:5,7]);
%         end
        switch GenTracker{j}
            case 'KO'
                KOfrac{i} = [KOfrac{i}; [active_frac,suppress_frac]];
            case 'HET'
                HETfrac{i} = [HETfrac{i}; [active_frac,suppress_frac]];
            case 'WT'
                WTfrac{i} = [WTfrac{i}; [active_frac,suppress_frac]];
            case 'C57'
                C57frac{i} = [C57frac{i};[active_frac,suppress_frac]];
        end
    end
    nbehav = size(C57frac{i},2);
    
    if i < 3
        xlabels = [strcat(toyBvs,'-act'),strcat(toyBvs,'-sup')];
    else
        xlabels = [strcat(catBvs,'-act'),strcat(catBvs,'-sup')];
    end

    XZBoxPlot(C57frac(i),[ones(1,nbehav/2),ones(1,nbehav/2)*2],[1:2:nbehav, 2:2:nbehav],xlabels,[]);
    xtickangle(45); set(gca,'fontsize',9,'TickLength',[0.004,0.004]); title(types{i}); legend({'suppress','active'});
%     figure, a = subplot(3,1,1);
%     XZBoxPlot(KOfrac(i),[ones(1,nbehav/2),ones(1,nbehav/2)*2],[1:2:nbehav, 2:2:nbehav],xlabels,[],a);
%     xtickangle(45);set(a,'fontsize',9,'TickLength',[0.004,0.004]);
%     title('KO'); legend({'suppress','active'});
%     a = subplot(3,1,2);
%     XZBoxPlot(HETfrac(i),[ones(1,nbehav/2),ones(1,nbehav/2)*2],[1:2:nbehav, 2:2:nbehav],xlabels,[],a);
%     xtickangle(45); set(a,'fontsize',9,'TickLength',[0.004,0.004]);
%     title('HET');legend({'suppress','active'});
%     a = subplot(3,1,3);
%     XZBoxPlot(WTfrac(i),[ones(1,nbehav/2),ones(1,nbehav/2)*2],[1:2:nbehav, 2:2:nbehav],xlabels,[],a);
%     xtickangle(45);set(a,'fontsize',9); set(a,'fontsize',9,'TickLength',[0.004,0.004]);
%     title('WT');legend({'suppress','active'});
%     sgtitle(types{i});
end


%% Plot fraction of cells of ROC output
types = {'self','partner','selfdrvd','partnerdrvd'};
gens = {'KO','HET','WT'};
PLSCthreshold = 0.95;
 % first half active second half suppress
for i = 1:numel(types)
    KOactive = []; HETactive = []; WTactive = [];
    for j = 1:numel(all_active.(types{i}))
        active_map1 = all_active.(types{i}){j};
% get PLSC significant cells ID

        M1_trace = allPairs{PairTracker(ceil(j/2))}{1}; M2_trace = allPairs{PairTracker(ceil(j/2))}{2};
        sessionN = find(~(cellfun(@isempty,M1_trace.Behavior)));
        M1_trace = zscore(M1_trace.MS{sessionN}.FiltTraces(M1_trace.TimeStamp.mapTs{sessionN}.M1toM2,:));
        M2_trace = zscore(M2_trace.MS{sessionN}.FiltTraces);
        [~,~,USV] = GetPLSC(M1_trace, M2_trace);
        M2_idx = find(or(USV.V(:,1)>=quantile(USV.V(:,1),PLSCthreshold), USV.V(:,1)<=quantile(USV.V(:,1),1-PLSCthreshold)));
        M1_idx = find(or(USV.U(:,1)>=quantile(USV.U(:,1),PLSCthreshold), USV.U(:,1)<=quantile(USV.U(:,1),1-PLSCthreshold)));

        if mod(j,2) % mouse 1
            active_map1 = active_map1(M1_idx,:);
        else
            active_map1 = active_map1(M2_idx,:);
        end

        ncell = size(active_map1,1);
        suppress_map1 = active_map1;
        active_map1(active_map1 < 0) = 0;
        suppress_map1(suppress_map1>0) = 0;
        active_map1 = round(sum(active_map1,1)/ncell,2);
        suppress_map1 = round(abs(sum(suppress_map1,1)/ncell),2);
        switch GenTracker{j}
            case 'KO'
                KOactive = [KOactive; [active_map1,suppress_map1,j]];
            case 'HET'
                HETactive = [HETactive; [active_map1, suppress_map1,j]];
            case 'WT'
                WTactive = [WTactive; [active_map1, suppress_map1,j]];
        end
    end
    if i < 3
        behavNames = rocStruct{1}.ROC.exp.other.self.Behavior.EventNames;
    else
        behavNames = rocStruct{1}.ROC.exp.other.derived_self.Behavior.EventNames;
    end
    prefix = 'PLSC1-';
    ncol = size(KOactive,2)-1;
    heatMap(KOactive(:,1:ncol/2),'xTickLabels',behavNames,'yTickLabels',ATracker(KOactive(:,end)),'fontSize',12);
    set(gcf,'Position',[0,0,1600,1000]); saveas(gcf,[targetFolder,'\',prefix,'KOactive-',types{i},'.png'],'png');
    heatMap(KOactive(:,ncol/2+1:ncol),'xTickLabels',behavNames,'yTickLabels',ATracker(KOactive(:,end)),'fontSize',12);
    set(gcf,'Position',[0,0,1600,1000]); saveas(gcf,[targetFolder,'\',prefix,'KOsuppress-',types{i},'.png'],'png');
    heatMap(HETactive(:,1:ncol/2),'xTickLabels',behavNames,'yTickLabels',ATracker(HETactive(:,end)),'fontSize',12);
    set(gcf,'Position',[0,0,1600,1000]); saveas(gcf,[targetFolder,'\',prefix,'HETactive-',types{i},'.png'],'png');
    heatMap(HETactive(:,ncol/2+1:ncol),'xTickLabels',behavNames,'yTickLabels',ATracker(HETactive(:,end)),'fontSize',12);
    set(gcf,'Position',[0,0,1600,1000]); saveas(gcf,[targetFolder,'\',prefix,'HETsuppress-',types{i},'.png'],'png');
    heatMap(WTactive(:,1:ncol/2),'xTickLabels',behavNames,'yTickLabels',ATracker(WTactive(:,end)),'fontSize',12);
    set(gcf,'Position',[0,0,1600,1000]); saveas(gcf,[targetFolder,'\',prefix,'WTactive-',types{i},'.png'],'png');
    heatMap(WTactive(:,ncol/2+1:ncol),'xTickLabels',behavNames,'yTickLabels',ATracker(WTactive(:,end)),'fontSize',12);
    set(gcf,'Position',[0,0,1600,1000]); saveas(gcf,[targetFolder,'\',prefix,'WTsuppress-',types{i},'.png'],'png');
    close all;
end

%% Plot by pairs
%% Plot fraction of cells of ROC output
types = {'self','partner','selfdrvd','partnerdrvd'};
gens = {'KO','HET','WT'};
PLSCthreshold = 0.0;
 % first half active second half suppress
for j = 1:2:numel(all_active.(types{1}))
    indi_bv = []; cat_bv = [];
    figure;
    for i = 1:numel(types)
        active_map1 = all_active.(types{i}){j};
        active_map2 = all_active.(types{i}){j+1};
% get PLSC significant cells ID
        M1 = allPairs{PairTracker(ceil(j/2))}{1}; M2 = allPairs{PairTracker(ceil(j/2))}{2};
        
        sessionN = find(~(cellfun(@isempty,M1.Behavior)));
        M1_trace = zscore(M1.MS{sessionN}.FiltTraces(M1.TimeStamp.mapTs{sessionN}.M1toM2,:));
        M2_trace = zscore(M2.MS{sessionN}.FiltTraces);
        [~,~,USV] = GetPLSC(M1_trace, M2_trace);
        M2_idx = find(or(USV.V(:,1)>=quantile(USV.V(:,1),PLSCthreshold), USV.V(:,1)<=quantile(USV.V(:,1),1-PLSCthreshold)));
        M1_idx = find(or(USV.U(:,1)>=quantile(USV.U(:,1),PLSCthreshold), USV.U(:,1)<=quantile(USV.U(:,1),1-PLSCthreshold)));

        
            active_map1 = active_map1(M1_idx,:);
        
            active_map2 = active_map2(M2_idx,:);
        

        ncell = size(active_map1,1);
        suppress_map1 = active_map1;
        active_map1(active_map1 < 0) = 0;
        suppress_map1(suppress_map1>0) = 0;
        active_map1 = round(sum(active_map1,1)/ncell,2);
        suppress_map1 = round(abs(sum(suppress_map1,1)/ncell),2);
        ncell = size(active_map2,1);
        suppress_map2 = active_map2;
        active_map2(active_map2 < 0) = 0;
        suppress_map2(suppress_map2>0) = 0;
        active_map2 = round(sum(active_map2,1)/ncell,2);
        suppress_map2 = round(abs(sum(suppress_map2,1)/ncell),2);
        if i <= 2
            indi_bv = [indi_bv;active_map1;active_map2;suppress_map1;suppress_map2];
        else
            cat_bv = [cat_bv; active_map1;active_map2;suppress_map1;suppress_map2];
        end
    end
    YtickLabel = {[M1.AnimalID,'-selfact'],[M2.AnimalID,'-selfact'],...
        [M1.AnimalID,'-selfsupp'],[M2.AnimalID,'-selfsupp'],...
        [M1.AnimalID,'-partnact'],[M2.AnimalID,'-partnact'],...
        [M1.AnimalID,'-partnsupp'],[M2.AnimalID,'-partnsupp']};
% prepare windowed behavior trace
    M1_B = M1.Behavior{sessionN}; M2_B = M2.Behavior{sessionN};
% determine the size of plot
    b_idx = find(~and(cellfun(@isempty,M1_B.OnsetTimes), cellfun(@isempty,M2_B.OnsetTimes)));
    total_b = numel(b_idx);
    total_panel = total_b + mod(total_b,2);
    
    subplot(3,total_panel,[1:round(total_panel*0.7)]);
    
    ExpID = strsplit(M1.ExperimentID,'_'); ExpID = ExpID{1};
    heatMap(indi_bv,'xTickLabels',rocStruct{1}.ROC.exp.other.self.Behavior.EventNames,'yTickLabels',YtickLabel,'fontSize',10,'colorbar','off');
    sgtitle([ExpID,'-',M1.AnimalID,M1.GenType,'-',M2.AnimalID,M2.GenType]);
%     axp = struct(gca); axp.Axes.XAxisLocation = 'top';
    subplot(3,total_panel,[round(total_panel*0.7)+1:total_panel]);
    heatMap(cat_bv,'xTickLabels',{'sniff','aggressive','defensive','exploratory','initiative','other'},'fontSize',10,'colorbar','off');
%     axp = struct(gca); axp.Axes.XAxisLocation = 'top';
    for k = 1:total_b
        M1_logi =M1_B.LogicalVecs{b_idx(k)}; M2_logi = M2_B.LogicalVecs{b_idx(k)};
        % get smooth trace
        M1_btrace = SmoothTrace(M1_logi,20*30,30);
        M2_btrace = SmoothTrace(M2_logi,20*30,30);
        subplot(3,total_panel,[2*(k-1)+total_panel+1,2*k+total_panel], 'NextPlot','add'); plot(M1_btrace,'r-'); plot(M2_btrace,'b-'); legend({M1.AnimalID, M2.AnimalID});
        title(M1_B.EventNames{b_idx(k)});
    end
    set(gcf,"Position",[-1919,41,1920,963]);
    saveas(gcf,[targetFolder,'ROC_byPair\',ExpID,'-',M1.AnimalID,M1.GenType,'-',M2.AnimalID,M2.GenType,'.png'],'png');
    close all;
    


    
   
end
%%
toyPexpP = cell(1,2); toyPexpN = cell(1,2); toyNexpP = cell(1,2); toyNexpN = cell(1,2); % first individual second derived
[all_active.covmat_indiv_pos, all_active.covmat_indiv_neg,all_active.covmat_indiv_pos1neg2,all_active.covmat_indiv_pos2neg1] = deal({});
[all_active.covmat_drvd_pos, all_active.covmat_drvd_neg,all_active.covmat_drvd_pos1neg2,all_active.covmat_drvd_pos2neg1] = deal({});
for i = 1 : 1/2*numel(all_active.self)
    temp_toy = all_active.self{i}; temp_exp = all_active.self{i+10};
    temp_toy_d = all_active.selfdrvd{i}; temp_exp_d = all_active.selfdrvd{i+10};
    [temp_toy_pos,temp_toy_neg] = deal(temp_toy); [temp_exp_pos,temp_exp_neg] = deal(temp_exp); 
    [temp_toy_d_pos,temp_toy_d_neg] = deal(temp_toy_d); [temp_exp_d_pos, temp_exp_d_neg] = deal(temp_exp_d);
    temp_toy_pos(temp_toy_pos < 0) = 0; temp_exp_pos(temp_exp_pos<0) = 0; 
    temp_exp_neg(temp_exp_neg > 0) = 0; temp_toy_neg(temp_toy_neg > 0) = 0;
    temp_toy_d_pos(temp_toy_d_pos < 0) = 0; temp_exp_d_pos(temp_exp_d_pos<0) = 0; 
    temp_exp_d_neg(temp_exp_d_neg > 0) = 0; temp_toy_d_neg(temp_toy_d_neg > 0) = 0;


    all_active.covmat_indiv_pos = [all_active.covmat_indiv_pos, temp_toy_pos' * temp_exp_pos];
    all_active.covmat_indiv_neg = [all_active.covmat_indiv_neg, temp_toy_neg' * temp_exp_neg];
    all_active.covmat_indiv_pos1neg2 = [all_active.covmat_indiv_pos1neg2, temp_toy_pos' * abs(temp_exp_neg)];
    all_active.covmat_indiv_pos2neg1 = [all_active.covmat_indiv_pos2neg1, abs(temp_toy_neg') * temp_exp_pos];
    figure; subplot(2,2,1); 
%     heatMap(temp_toy_pos' * temp_exp_pos,'xTickLabels',expBvs,'yTickLabels',toyBvs,'fontSize',10); title('toyPos-expPos'); subplot(2,2,2); 
%     heatMap(temp_toy_neg' * temp_exp_neg,'xTickLabels',expBvs,'yTickLabels',toyBvs,'fontSize',10); title('toyNeg-expNeg'); subplot(2,2,3);
%     heatMap(temp_toy_pos' * abs(temp_exp_neg),'xTickLabels',expBvs,'yTickLabels',toyBvs,'fontSize',10); title('toyPos-expNeg'); subplot(2,2,4);
%     heatMap(abs(temp_toy_neg') * temp_exp_pos,'xTickLabels',expBvs,'yTickLabels',toyBvs,'fontSize',10); title('toyNeg-expPos');
    toyPexpP{1} = cat(3,toyPexpP{1},temp_toy_pos' * temp_exp_pos/size(temp_exp_pos,1));
    toyNexpP{1} = cat(3,toyNexpP{1},abs(temp_toy_neg)' * temp_exp_pos/size(temp_exp_pos,1));
    toyPexpN{1} = cat(3,toyPexpN{1},temp_toy_pos' * abs(temp_exp_neg)/size(temp_exp_pos,1));
    toyNexpN{1} = cat(3,toyNexpN{1},temp_toy_neg' * temp_exp_neg/size(temp_exp_pos,1));
    
    PairID = strsplit(ATracker{2*i-1},'-'); PairID = PairID{1};
%     set(gcf,'Position',[0,0,1920,1080]); saveas(gcf,[targetFolder,'\',PairID,'1.png'],'png');
    close all;
    all_active.covmat_drvd_pos = [all_active.covmat_drvd_pos, temp_toy_d_pos' * temp_exp_d_pos];
    all_active.covmat_drvd_neg = [all_active.covmat_drvd_neg, temp_toy_d_neg' * temp_exp_d_neg];
    all_active.covmat_drvd_pos1neg2 = [all_active.covmat_drvd_pos1neg2, temp_toy_d_pos' * abs(temp_exp_d_neg)];
    all_active.covmat_drvd_pos2neg1 = [all_active.covmat_drvd_pos2neg1, abs(temp_toy_d_neg') * temp_exp_d_pos];
    figure; subplot(2,2,1); 
%     heatMap(temp_toy_d_pos' * temp_exp_d_pos,'xTickLabels',expdrvBvs,'yTickLabels',toydrvBvs,'fontSize',10); title('toyPos-expPos'); subplot(2,2,2); 
%     heatMap(temp_toy_d_neg' * temp_exp_d_neg,'xTickLabels',expdrvBvs,'yTickLabels',toydrvBvs,'fontSize',10); title('toyNeg-expNeg'); subplot(2,2,3);
%     heatMap(temp_toy_d_pos' * abs(temp_exp_d_neg),'xTickLabels',expdrvBvs,'yTickLabels',toydrvBvs,'fontSize',10); title('toyPos-expNeg'); subplot(2,2,4);
%     heatMap(abs(temp_toy_d_neg') * temp_exp_d_pos,'xTickLabels',expdrvBvs,'yTickLabels',toydrvBvs,'fontSize',10); title('toyNeg-expPos');
    toyPexpP{2} = cat(3,toyPexpP{2},temp_toy_d_pos' * temp_exp_d_pos/size(temp_exp_d_pos,1));
    toyNexpP{2} = cat(3,toyNexpP{2},abs(temp_toy_d_neg)' * temp_exp_d_pos/size(temp_exp_d_pos,1));
    toyPexpN{2} = cat(3,toyPexpN{2},temp_toy_d_pos' * abs(temp_exp_d_neg)/size(temp_exp_d_pos,1));
    toyNexpN{2} = cat(3,toyNexpN{2},temp_toy_d_neg' * temp_exp_d_neg/size(temp_exp_d_pos,1));
    PairID = strsplit(ATracker{2*i-1},'-'); PairID = PairID{1};
%     set(gcf,'Position',[0,0,1920,1080]); saveas(gcf,[targetFolder,'\',PairID,'2.png'],'png');
    close all
    
end

figure; subplot(2,2,1); 
    heatMap(round(mean(toyPexpP{1}(1:end-1,1:end-1,:),3,'omitnan'),2),'xTickLabels',expBvs(1:end-1),'yTickLabels',toyBvs(1:end-1),'fontSize',10); title('toyPos-expPos'); subplot(2,2,2); 
    heatMap(round(mean(toyNexpN{1}(1:end-1,1:end-1,:),3,'omitnan'),2),'xTickLabels',expBvs(1:end-1),'yTickLabels',toyBvs(1:end-1),'fontSize',10); title('toyNeg-expNeg'); subplot(2,2,3);
    heatMap(round(mean(toyPexpN{1}(1:end-1,1:end-1,:),3,'omitnan'),2),'xTickLabels',expBvs(1:end-1),'yTickLabels',toyBvs(1:end-1),'fontSize',10); title('toyPos-expNeg'); subplot(2,2,4);
    heatMap(round(mean(toyNexpP{1}(1:end-1,1:end-1,:),3,'omitnan'),2),'xTickLabels',expBvs(1:end-1),'yTickLabels',toyBvs(1:end-1),'fontSize',10); title('toyNeg-expPos');


figure; subplot(2,2,1); 
    heatMap(round(mean(toyPexpP{2}(3:6,1:5,:),3,'omitnan'),2),'xTickLabels',catBvs(1:5),'yTickLabels',catBvs(3:6),'fontSize',10); title('toyPos-expPos'); subplot(2,2,2); 
    heatMap(round(mean(toyNexpN{2}(3:6,1:5,:),3,'omitnan'),2),'xTickLabels',catBvs(1:5),'yTickLabels',catBvs(3:6),'fontSize',10); title('toyNeg-expNeg'); subplot(2,2,3);
    heatMap(round(mean(toyPexpN{2}(3:6,1:5,:),3,'omitnan'),2),'xTickLabels',catBvs(1:5),'yTickLabels',catBvs(3:6),'fontSize',10); title('toyPos-expNeg'); subplot(2,2,4);
    heatMap(round(mean(toyNexpP{2}(3:6,1:5,:),3,'omitnan'),2),'xTickLabels',catBvs(1:5),'yTickLabels',catBvs(3:6),'fontSize',10); title('toyNeg-expPos');

