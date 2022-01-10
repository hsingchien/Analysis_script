% find ROC cell identity
% A1
roc_test = rocStruct{1, 1}{1, 2}.ROC{1, 2}.all.derived_partner.Encoding;
roc = roc_test.auROC;
rocrand = roc_test.auROCrandDist;
nonNA_idx = find(~isnan(roc(1,:)));
nonNA_behav = roc_test.EventNames(nonNA_idx);
higher_bounds = ones(size(roc));
lower_bounds = zeros(size(roc));
for i = nonNA_idx
   for j = 1:size(roc,1)
       higher_bound = quantile(rocrand{j,i},0.975);
       lower_bound = quantile(rocrand{j,i},0.025);
       higher_bounds(j,i) = higher_bound;
       lower_bounds(j,i) = lower_bound;
   end
    
end

poscellID = {};
negcellID = {};
for i = nonNA_idx
    poscellID = [poscellID, find(roc(:,i)>=higher_bounds(:,i))];
    negcellID = [negcellID, find(roc(:,i)<=lower_bounds(:,i))];
end
