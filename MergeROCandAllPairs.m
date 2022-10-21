%% load
allExpID = {};
for k = 1:numel(allPairs)
    allExpID = [allExpID,allPairs{k}{1}.ExperimentID];
end
for i = 1:numel(flist)
    load(flist{i});
    fprintf('%s\n',flist{i});
    M1 = rocStruct{1}; M2 = rocStruct{2};
    EName = M1.ExperimentID; EName = strsplit(EName,'_');EName = [EName{1},'_'];
    Eid = find(contains(allExpID,EName));
    if isempty(M1.ROC)
        continue
    end
    thisPair = allPairs{Eid};
    aIDs = {thisPair{1}.AnimalID,thisPair{2}.AnimalID};
    A1 = find(strcmp(aIDs,M1.AnimalID)); A2 = find(strcmp(aIDs, M2.AnimalID));

    sessionID = find(~cellfun(@isempty,thisPair{1}.Behavior));

    M1 = M1.ROC.exp.other; M2 = M2.ROC.exp.other;
    f = fields(M1);
    for ff = 1:numel(f)
        if size(thisPair{A1}.MS{sessionID}.FiltTraces,1) == size(M1.(f{ff}).zCalciumTraces,1)
            M1.(f{ff}) = rmfield(M1.(f{ff}),{'zCalciumTraces','Behavior'});
            M2.(f{ff}) = rmfield(M2.(f{ff}),{'zCalciumTraces','Behavior'});
            thisPair{A1}.Behavior{sessionID}.ROC = M1;
            thisPair{A2}.Behavior{sessionID}.ROC = M2;
            fprintf('success\n');
        else
            fprintf('failed on %s\n', flist{i});
        end
    end

    allPairs{Eid} = thisPair;
    
    

end
