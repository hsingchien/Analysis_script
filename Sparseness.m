%% sparseness
% percentage of off cells

cur_cell = DLX_shank;
cur_type = 'DLX_shank';
social_bv = {'attack','chasing','escape','defend','flinch','tussling','threaten',...
    'general-sniffing','sniff_face','sniff_genital',...
    'approach','follow','interaction','attention','socialgrooming'}; 
sample_size = 50;
nonsocial_bv = {'dig','climb','exploreobj','stand'};
idle_bv = {'selfgrooming','other'};

for i = 1:length(cur_cell)
    for j = 1:2
        cur_animal = cur_cell{i}{j};
        for m = 1:length(cur_animal.MS)
            if isempty(cur_animal.Behavior{m})
                continue;
            end
%             cur_S = transpose(cur_animal.MS{m}.S);
            cur_S = cur_animal.MS{m}.S;
%             cur_S = cur_S(:,logical(cur_animal.MS{m}.goodCellVec));
            % random sample 50
%             cur_S = datasample(cur_S, sample_size, 2);
            cur_S_threshold = mean(cur_S,1) + 2*std(cur_S);
            cur_S_filtered = cur_S;
            cur_S_binary = cur_S >= repmat(cur_S_threshold,[size(cur_S,1),1]);
            cur_S_filtered(~cur_S_binary) = 0;
            cur_S_binary = double(cur_S_binary);
            twindow = 30;
            step = round(twindow/2);
            sparsenesses = [];
            pcts = [];
            tws = [];
            B = cur_animal.Behavior{m};
            Bt = cur_animal.TimeStamp.Ts{m}.Bv;
            Mst = cur_animal.TimeStamp.Ts{m}.Ms;
            for tw = 1:step:(size(cur_S_filtered,1)-twindow+1)
                cur_tw = cur_S_filtered(tw:(tw+twindow-1),:);
                cur_tw_bin = cur_S_binary(tw:(tw+twindow-1),:);
                pct = sum((sum(cur_tw_bin,1)>0))/size(cur_S_filtered,2);
                fr = sum(cur_tw,1); % 1 x n_neuron
                sparseness = (sum(fr)/size(cur_S_filtered,2))^2/sum(fr.^2/size(cur_S_filtered,2));
                tws = [tws,tw];
                pcts = [pcts, pct];
                sparsenesses = [sparsenesses, sparseness];
            end
            tws = tws + step;
            f = figure;
            a1 = subplot(2,1,1);a1.NextPlot = 'add';
            a2 = subplot(2,1,2);a2.NextPlot = 'add';
            plot(a1, Mst(tws)/1000, pcts); title(a1,'PCT');
            a1.XLim = [0, Mst(end)/1000];
            
            plot(a2, Mst(tws)/1000, sparsenesses); title(a2,'Sparseness');
            a2.XLim = [0, Mst(end)/1000];
            if ~isempty(B)
%                 AddBehavBlocks(a2,B,Bt/1000,social_bv,assignColors(social_bv));
                AddBehavBlocks(a2,B,Bt/1000,[nonsocial_bv,social_bv],[repmat([0,1,1],[length(nonsocial_bv),1]);repmat([1,0.6,0.6],[length(social_bv),1])]);
            end
            savefig(f,[cur_type,'\P',num2str(i),'A',num2str(j),cur_animal.videoInfo.session{m},'_',num2str(m)]);
            close(f);
        end
    end
end
%% sparseness average over behaviors
cur_cell = DLX_shank;
cur_type = 'DLX_shank';
nonsocial_sparse = [];
social_sparse = [];
other_sparse = [];
GenTracker = {};
% sample_size = 50;
for i = 1:length(cur_cell)
    for j = 1:2
        cur_animal = cur_cell{i}{j};
        for m = 1:length(cur_animal.MS)
            if isempty(cur_animal.Behavior{m})
                continue;
            end
%             for kk = 1:50
%                 cur_S = transpose(cur_animal.MS{m}.S);
                cur_S = cur_animal.MS{m}.S;
%                 cur_S = cur_S(:,logical(cur_animal.MS{m}.goodCellVec));
%                 cur_S = datasample(cur_S, sample_size, 2);
                cur_S_threshold = mean(cur_S,1) + 2*std(cur_S);
                cur_S_filtered = cur_S;
                cur_S_binary = cur_S >= repmat(cur_S_threshold,[size(cur_S,1),1]);
                cur_S_filtered(~cur_S_binary) = 0;
                cur_S_binary = double(cur_S_binary);
                twindow = 30;
                step = round(twindow/2);
                sparsenesses = [];
                pcts = [];
                tws = [];
                B = cur_animal.Behavior{m};
                Bt = cur_animal.TimeStamp.Ts{m}.Bv;
                Mst = cur_animal.TimeStamp.Ts{m}.Ms;
                for tw = 1:step:(size(cur_S_filtered,1)-twindow+1)
                    cur_tw = cur_S_filtered(tw:(tw+twindow-1),:);
                    cur_tw_bin = cur_S_binary(tw:(tw+twindow-1),:);
                    pct = sum((sum(cur_tw_bin,1)>0))/size(cur_S_filtered,2);
                    fr = sum(cur_tw,1); % 1 x n_neuron
                    sparseness = (sum(fr)/size(cur_S_filtered,2))^2/sum(fr.^2/size(cur_S_filtered,2));
                    tws = [tws,tw];
                    pcts = [pcts, pct];
                    sparsenesses = [sparsenesses, sparseness];
                end
                tws = tws + round(step/2); % take the middle
                behav_vecs = zeros([1,length(cur_animal.Behavior{m}.LogicalVecs{1})]);
                for ii = 1:length(cur_animal.Behavior{m}.LogicalVecs)
                    behav_vecs = behav_vecs + ii*cur_animal.Behavior{m}.LogicalVecs{ii};
                end
                behav_vecs = behav_vecs(cur_animal.TimeStamp.mapTs{m}.B2M);
                behav_vecs = behav_vecs(tws);
                [~,social_index] = ismember(social_bv, cur_animal.Behavior{m}.EventNames);
                [~,nonsocial_index] = ismember(nonsocial_bv, cur_animal.Behavior{m}.EventNames);
                [~,other_index] = ismember({'other'}, cur_animal.Behavior{m}.EventNames);
                social_tlogic = ismember(behav_vecs, social_index);
                nonsocial_tlogic = ismember(behav_vecs, nonsocial_index);
                other_tlogic = ismember(behav_vecs, other_index);
                temp1 = sparsenesses(social_tlogic);
                social_sparse = [social_sparse,mean(temp1(~isnan(temp1)))];
                temp2 = sparsenesses(nonsocial_tlogic);
                nonsocial_sparse = [nonsocial_sparse, mean(temp2(~isnan(temp2)))];
                temp3 = sparsenesses(other_tlogic);
                other_sparse = [other_sparse, mean(temp3(~isnan(temp3)))];
                GenTracker = [GenTracker, cur_cell{i}{j}.GenType];
%             end
        end
    end
end
% boxplot([social_sparse,nonsocial_sparse,other_sparse],[ones(1,length(social_sparse)),ones(1,length(nonsocial_sparse))*2, ones(1,length(other_sparse))*3])
other_sparse = other_sparse'; social_sparse = social_sparse'; nonsocial_sparse = nonsocial_sparse';
figure, a = subplot(1,3,1);
XZBoxPlot({other_sparse(strcmp(GenTracker,'KO')),other_sparse(strcmp(GenTracker,'HET')),other_sparse(strcmp(GenTracker,'WT'))}, ...
    1:3,1:3,{'KO','HET','WT'},[],a,true); title('Idle Sparseness'); ylim([0,0.24]); ylabel('Sparsenss (the lower the sparser)')
a = subplot(1,3,2);
XZBoxPlot({social_sparse(strcmp(GenTracker,'KO')),social_sparse(strcmp(GenTracker,'HET')),social_sparse(strcmp(GenTracker,'WT'))}, ...
    1:3,1:3,{'KO','HET','WT'},[],a,true); title('Social Sparseness');ylim([0,0.24]);
a = subplot(1,3,3);
XZBoxPlot({nonsocial_sparse(strcmp(GenTracker,'KO')),nonsocial_sparse(strcmp(GenTracker,'HET')),nonsocial_sparse(strcmp(GenTracker,'WT'))}, ...
    1:3,1:3,{'KO','HET','WT'},[],a,true); title('Nonsocial Sparseness');ylim([0,0.24]);

XZBoxPlot({other_sparse(strcmp(GenTracker,'KO')),other_sparse(strcmp(GenTracker,'HET')),other_sparse(strcmp(GenTracker,'WT')), ...
    nonsocial_sparse(strcmp(GenTracker,'KO')),nonsocial_sparse(strcmp(GenTracker,'HET')),nonsocial_sparse(strcmp(GenTracker,'WT')),...
    social_sparse(strcmp(GenTracker,'KO')),social_sparse(strcmp(GenTracker,'HET')),social_sparse(strcmp(GenTracker,'WT'))}, ...
    [1,2,3,1,2,3,1,2,3],[1,4,7,2,5,8,3,6,9],{'KO-idle','HET-idle','WT-idle','KO-nonsocial','HET-nonsocial','WT-nonsocial','KO-social','HET-social','WT-social'}, ...
    [],[],true); title('Sparseness'); ylim([0,0.24]); ylabel('Sparsenss (the lower the sparser)');
%%
cur_cell = CMK;
cur_type = 'CMK';

for i = 1:length(cur_cell)
    for j = 1:2
        cur_animal = cur_cell{i}{j};
        pairwise_toy = [];
        pairwise_exp = [];
        for m = 1:length(cur_animal.MS)
            if contains(cur_animal.videoInfo.session{m},'sep');
                continue;
            end
            if size(cur_animal.MS{m}.FiltTraces,1) < 8000
                continue;
            end
            covmat = reshape(corr(zscore(cur_animal.MS{m}.FiltTraces(:, cur_animal.MS{m}.goodCellVec))),[],1);
            if contains(cur_animal.videoInfo.session{m},'toy') & isempty(pairwise_toy)
                pairwise_toy = covmat;
            elseif contains(cur_animal.videoInfo.session{m},'exp') & isempty(pairwise_exp)
                pairwise_exp = covmat;
            end
        end
        
        figure, plot(pairwise_toy(pairwise_toy~=1), pairwise_exp(pairwise_toy~=1),'b.');
    end
end




                