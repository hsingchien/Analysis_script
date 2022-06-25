% get pca from toy/sep session and project exp on it, count the variance
% explained
% mean bin trace to 1s. 
% factors: exp,toy,sep; gentype: ko,wt,het;
gens = {'KO','HET','WT'};
sep_var = {[],[],[]}; % KO, HET, WT
toy_var = {[],[],[]};
exp_sep_var = {[],[],[]};
exp_toy_var = {[],[],[]};
exp_sep_cos = {[],[],[]};
exp_sep_pls_cos = exp_sep_cos;
exp_toy_cos = {[],[],[]};
exp_toy_pls_cos = exp_toy_cos;
Aid = {[],[],[]};
pcn = 1;
for i = 1:length(allPairs)
    m1 = allPairs{i}{1};
    m2 = allPairs{i}{2};
    if strcmp(m1.GenType,'C57') | strcmp(m2.GenType,'C57')
        continue
    end
    % pca sep
    sepsessions = find(contains(m1.videoInfo.session,'sep'));
    expsessions = find(contains(m1.videoInfo.session,'exp'));
    toysessions = find(contains(m1.videoInfo.session,'toy'));
    if isempty(toysessions)
        continue
    end
    sepsessions = sepsessions(end); % only keep the last one
    Ftrace1_sep = zscore(BinData(m1.MS{sepsessions}.FiltTraces,15,'mean'));
    Ftrace2_sep = zscore(BinData(m2.MS{sepsessions}.FiltTraces,15,'mean'));
    [sep_coeff1,~,~,~,sep_explained1] = pca(Ftrace1_sep);
    [sep_coeff2,~,~,~,sep_explained2] = pca(Ftrace2_sep);
    [~,~,USV_sep] = GetPLSC(Ftrace1_sep, Ftrace2_sep);
    exp_explained1 = 0;
    exp_explained2 = 0;
    exp_cos1 = 0; exp_cos2 = 0; sep_cos1 = 0; sep_cos2 = 0;
    totall1 = 0; totall2 = 0;
    for j = 1:numel(expsessions)
        Ftrace1_exp = zscore(BinData(m1.MS{expsessions(j)}.FiltTraces(m1.TimeStamp.mapTs{expsessions(j)}.M1toM2,:),15,'mean'));
        Ftrace2_exp = zscore(BinData(m2.MS{expsessions(j)}.FiltTraces,15,'mean'));
        Ftrace1_proj = Ftrace1_exp*sep_coeff1;
        Ftrace2_proj = Ftrace2_exp*sep_coeff2;
        var_1 = diag(cov(Ftrace1_proj));
        var_2 = diag(cov(Ftrace2_proj));
        exp_explained1 = exp_explained1 + sum(var_1(pcn))/sum(var_1)*size(Ftrace1_exp,1); totall1 = totall1+size(Ftrace1_exp,1);
        exp_explained2 = exp_explained2 + sum(var_2(pcn))/sum(var_2)*size(Ftrace2_exp,1); totall2 = totall2+size(Ftrace2_exp,1);
        [~,~,USV] = GetPLSC(Ftrace1_exp,Ftrace2_exp);
        exp_cos1 = exp_cos1 + dot(USV.U(:,1),sep_coeff1(:,1))/(norm(USV.U(:,1))*norm(sep_coeff1(:,1)))*size(Ftrace1_exp,1);
        exp_cos2 = exp_cos2 + dot(USV.V(:,1),sep_coeff2(:,1))/(norm(USV.V(:,1))*norm(sep_coeff2(:,1)))*size(Ftrace2_exp,1);
        sep_cos1 = sep_cos1 + dot(USV.U(:,1), USV_sep.U(:,1))/(norm(USV.U(:,1))*norm(USV_sep.U(:,1)))*size(Ftrace1_exp,1);
        sep_cos2 = sep_cos2 + dot(USV.V(:,1), USV_sep.V(:,1))/(norm(USV.V(:,1))*norm(USV_sep.V(:,1)))*size(Ftrace2_exp,1);
    end
    exp_explained1 = exp_explained1/totall1;
    exp_explained2 = exp_explained2/totall2;
    sep_cos1 = sep_cos1/totall1; sep_cos2 = sep_cos2/totall2;
    exp_cos1 = exp_cos1/totall1; exp_cos2 = exp_cos2/totall2;
    exp_sep_var{find(strcmp(gens,m1.GenType))} = [exp_sep_var{find(strcmp(gens,m1.GenType))}; exp_explained1];
    exp_sep_var{find(strcmp(gens,m2.GenType))} = [exp_sep_var{find(strcmp(gens,m2.GenType))}; exp_explained2];
    sep_var{find(strcmp(gens,m1.GenType))} = [sep_var{find(strcmp(gens,m1.GenType))}; sum(sep_explained1(pcn))/100];
    sep_var{find(strcmp(gens,m2.GenType))} = [sep_var{find(strcmp(gens,m2.GenType))}; sum(sep_explained2(pcn))/100];
    exp_sep_cos{find(strcmp(gens,m1.GenType))} = [exp_sep_cos{find(strcmp(gens,m1.GenType))}; acos(abs(exp_cos1))];
    exp_sep_cos{find(strcmp(gens,m2.GenType))} = [exp_sep_cos{find(strcmp(gens,m2.GenType))}; acos(abs(exp_cos2))];
    exp_sep_pls_cos{find(strcmp(gens,m1.GenType))} = [exp_sep_pls_cos{find(strcmp(gens,m1.GenType))}; acos(abs(sep_cos1))];
    exp_sep_pls_cos{find(strcmp(gens,m2.GenType))} = [exp_sep_pls_cos{find(strcmp(gens,m2.GenType))}; acos(abs(sep_cos2))];

    
    Aid{find(strcmp(gens,m1.GenType))} = [Aid{find(strcmp(gens,m1.GenType))};10*i+1];
    Aid{find(strcmp(gens,m2.GenType))} = [Aid{find(strcmp(gens,m2.GenType))};10*i+2];
    if ~isempty(toysessions)
        toysessions = toysessions(end); % only keep the last one
        Ftrace1_toy = zscore(BinData(m1.MS{toysessions}.FiltTraces,15,'mean'));
        Ftrace2_toy = zscore(BinData(m2.MS{toysessions}.FiltTraces,15,'mean'));
        [toy_coeff1,~,~,~,toy_explained1] = pca(Ftrace1_toy);
        [toy_coeff2,~,~,~,toy_explained2] = pca(Ftrace2_toy);
        [~,~,USV_toy] = GetPLSC(Ftrace1_toy,Ftrace2_toy);
        exp_explained1 = 0;
        exp_explained2 = 0;
        totall1 = 0; totall2 = 0;
        exp_cos1 = 0; exp_cos2 = 0; toy_cos1 = 0; toy_cos2 = 0;
        for j = 1:numel(expsessions)
            Ftrace1_exp = zscore(BinData(m1.MS{expsessions(j)}.FiltTraces(m1.TimeStamp.mapTs{expsessions(j)}.M1toM2,:),15,'mean'));
            Ftrace2_exp = zscore(BinData(m2.MS{expsessions(j)}.FiltTraces,15,'mean'));
            Ftrace1_proj = Ftrace1_exp*toy_coeff1;
            Ftrace2_proj = Ftrace2_exp*toy_coeff2;
            var_1 = diag(cov(Ftrace1_proj));
            var_2 = diag(cov(Ftrace2_proj));
            exp_explained1 = exp_explained1 + sum(var_1(pcn))/sum(var_1)*size(Ftrace1_exp,1); totall1 = totall1+size(Ftrace1_exp,1);
            exp_explained2 = exp_explained2 + sum(var_2(pcn))/sum(var_2)*size(Ftrace2_exp,1); totall2 = totall2+size(Ftrace2_exp,1);
            [~,~,USV] = GetPLSC(Ftrace1_exp,Ftrace2_exp);
            exp_cos1 = exp_cos1 + dot(USV.U(:,1),toy_coeff1(:,1))/(norm(USV.U(:,1))*norm(toy_coeff1(:,1)))*size(Ftrace1_exp,1);
            exp_cos2 = exp_cos2 + dot(USV.V(:,1),toy_coeff2(:,1))/(norm(USV.V(:,1))*norm(toy_coeff2(:,1)))*size(Ftrace2_exp,1);
            toy_cos1 = toy_cos1 + dot(USV.U(:,1), USV_toy.U(:,1))/(norm(USV.U(:,1))*norm(USV_toy.U(:,1)))*size(Ftrace1_exp,1);
            toy_cos2 = toy_cos2 + dot(USV.V(:,1), USV_toy.V(:,1))/(norm(USV.V(:,1))*norm(USV_toy.V(:,1)))*size(Ftrace2_exp,1);
        end
        exp_cos1 = exp_cos1/totall1; exp_cos2 = exp_cos2/totall2;
        toy_cos1 = toy_cos1/totall1; toy_cos2 = toy_cos2/totall2;
        exp_explained1 = exp_explained1/totall1;
        exp_explained2 = exp_explained2/totall2;
        exp_toy_var{find(strcmp(gens,m1.GenType))} = [exp_toy_var{find(strcmp(gens,m1.GenType))}; exp_explained1];
        exp_toy_var{find(strcmp(gens,m2.GenType))} = [exp_toy_var{find(strcmp(gens,m2.GenType))}; exp_explained2];
        toy_var{find(strcmp(gens,m1.GenType))} = [toy_var{find(strcmp(gens,m1.GenType))}; sum(toy_explained1(pcn))/100];
        toy_var{find(strcmp(gens,m2.GenType))} = [toy_var{find(strcmp(gens,m2.GenType))}; sum(toy_explained2(pcn))/100];
        exp_toy_cos{find(strcmp(gens,m1.GenType))} = [exp_toy_cos{find(strcmp(gens,m1.GenType))}; acos(abs(exp_cos1))];
        exp_toy_cos{find(strcmp(gens,m2.GenType))} = [exp_toy_cos{find(strcmp(gens,m2.GenType))}; acos(abs(exp_cos2))];
        exp_toy_pls_cos{find(strcmp(gens,m1.GenType))} = [exp_toy_pls_cos{find(strcmp(gens,m1.GenType))}; acos(abs(toy_cos1))];
        exp_toy_pls_cos{find(strcmp(gens,m2.GenType))} = [exp_toy_pls_cos{find(strcmp(gens,m2.GenType))}; acos(abs(toy_cos2))];
        
    end

end
%%
% XZBoxPlot([sep_var,exp_sep_var],[1,2,3,1,2,3],[1,3,5,2,4,6],{'KO-sep','HET-sep','WT-sep','KO-exp','HET-exp','WT-exp'},[],1,true); ylim([0,0.4]);
% XZBoxPlot([toy_var,exp_toy_var],[1,2,3,1,2,3],[1,3,5,2,4,6],{'KO-sep','HET-sep','WT-sep','KO-exp','HET-exp','WT-exp'},[],1,true); ylim([0,0.4]);
XZBoxPlot([sep_var,exp_sep_var,exp_toy_var,toy_var], ...
    [1,2,3,1,2,3,1,2,3,1,2,3],[1,5,9,2,6,10,3,7,11,4,8,12], ...
    {'KO-sep','HET-sep','WT-sep','KO-exp-sep','HET-exp-sep','WT-exp-sep','KO-exp-toy','HET-exp-toy','WT-exp-toy','KO-toy','HET-toy','WT-toy'},[],1,true);
title('frac variance explained by PC1')
XZBoxPlot([exp_sep_cos,exp_toy_cos],[1,2,3,1,2,3],[1,3,5,2,4,6],{'KO-exp-sep','HET-exp-sep','WT-exp-sep','KO-exp-toy','HET-exp-toy','WT-exp-toy'},[],1,true);
title('Cos between PLSC1 and sep or toy PC1');
XZBoxPlot([exp_sep_pls_cos,exp_toy_pls_cos],[1,2,3,1,2,3],[1,3,5,2,4,6],{'KO-exp-sep','HET-exp-sep','WT-exp-sep','KO-exp-toy','HET-exp-toy','WT-exp-toy'},[],1,true);
title('Cos between exp PLSC1 and sep or toy PLSC1');