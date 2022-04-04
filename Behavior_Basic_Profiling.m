ds_behav = [1,2,3,5,6,7,9,10];
all_behav_exp = {'attack','chasing','tussling','threaten','escape','defend',...
    'flinch','general-sniffing','sniff_face','sniff_genital','approach',...
    'follow','interaction', 'socialgrooming', 'mount','dig',...
    'selfgrooming', 'climb', 'exploreobj', 'biteobj', 'stand', 'nesting','human_interfere', 'other'};
social_bv = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
nonsocial_bv = [16,17,18,19,20,21,22,23];
agonist_bv = [1,2,3,4];
defensive_bv = [5,6,7];
nondefensive_social_bv = setdiff(social_bv, defensive_bv);
sniff_bv = [8,9,10];
initiative_bv = [11,12];
allbehav_idx = {social_bv,nonsocial_bv,agonist_bv,defensive_bv,nondefensive_social_bv,sniff_bv,initiative_bv};
hm_interfere = 23;
%% percentage of each behavior category
categorical_bv = cell(1,length(allbehav_idx));
individual_bv = cell(1,22);
animalID = {};
animalGen = {};

for i = ds_behav
    M1 = allPairs{i}{1};
    M2 = allPairs{i}{2};
    Gen1 = M1.GenType;
    Gen2 = M2.GenType;
    % remove human_interfere
    non_human = ~(M1.Behavior{2}.LogicalVecs{hm_interfere} + M2.Behavior{2}.LogicalVecs{hm_interfere});
    animalGen = [animalGen;Gen1;Gen2];
    % deal with behavior categories
    for j = 1:length(allbehav_idx)
        cur_behav_idx = allbehav_idx{j};
        collapsed_logical_1 = zeros(1,sum(non_human));
        collapsed_logical_2 = collapsed_logical_1;
        for k = cur_behav_idx
            collapsed_logical_1 = collapsed_logical_1 + M1.Behavior{2}.LogicalVecs{k}(non_human);
            collapsed_logical_2 = collapsed_logical_2 + M2.Behavior{2}.LogicalVecs{k}(non_human);
        end
        perct1 = sum(collapsed_logical_1) / length(collapsed_logical_1);
        perct2 = sum(collapsed_logical_2) / length(collapsed_logical_2);
        categorical_bv{j} = [categorical_bv{j};perct1;perct2];
    end
    % deal with individual behaviors
    for j = 1:22 % except human & other
       collapsed_perct_1 = sum(M1.Behavior{2}.LogicalVecs{j}(non_human))/sum(non_human);
       collapsed_perct_2 = sum(M2.Behavior{2}.LogicalVecs{j}(non_human))/sum(non_human);
       individual_bv{j} = [individual_bv{j};collapsed_perct_1;collapsed_perct_2];
    end
    animalID = [animalID;M1.AnimalID;M2.AnimalID];
    aimalGen = [animalGen;M1.GenType;M2.GenType];
    
end
%% generate table, vars = 6 behav types + 1 gen type + aID
T1 = table(categorical_bv{:}, animalID, animalGen, 'VariableNames',{'social','nonsocial','agonist','defensive','nondefensive_social','sniff','initiative','animalID','animalGen'});
writetable(T1,'categorical_bv_prct.csv');
T2 = table(individual_bv{:}, animalID, animalGen, 'VariableNames',[all_behav_exp(1:22),'animalID','animalGen']);
writetable(T2,'individual_bv_prct.csv');

        
        