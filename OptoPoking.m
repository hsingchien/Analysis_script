[B1_1,~] = BehavStruExtract('04-14-48.txt',1);
[B2_1,~] = BehavStruExtract('04-14-48.txt',2);
[B1_2,~] = BehavStruExtract('04-35-49.txt',1);
[B2_2,~] = BehavStruExtract('04-35-49.txt',2);
[B1_3,~] = BehavStruExtract('04-44-58.txt',1);
[B2_3,~] = BehavStruExtract('04-44-58.txt',2);
[opto1_2,~] = BehavStruExtract('04-35-49.txt',3);
[opto2_2,~] = BehavStruExtract('04-35-49.txt',4);
[opto1_3,~] = BehavStruExtract('04-44-58.txt',3);
[opto2_3,~] = BehavStruExtract('04-44-58.txt',4);
%%
poken = 1; sniffn = 4; opton = 12;
B1_1_log = B1_1.LogicalVecs{poken}+B1_1.LogicalVecs{sniffn};
B2_1_log = B2_1.LogicalVecs{poken}+B2_1.LogicalVecs{sniffn};
B1_2_log = B1_2.LogicalVecs{poken}+B1_2.LogicalVecs{sniffn};
B2_2_log = B2_2.LogicalVecs{poken}+B2_2.LogicalVecs{sniffn};
B1_3_log = B1_3.LogicalVecs{poken}+B1_3.LogicalVecs{sniffn};
B2_3_log = B2_3.LogicalVecs{poken}+B2_3.LogicalVecs{sniffn};

B1_2_opto = [find(opto1_2.LogicalVecs{opton},1), find(opto1_2.LogicalVecs{opton},1,'last')];
B2_2_opto = [find(opto2_2.LogicalVecs{opton},1), find(opto2_2.LogicalVecs{opton},1,'last')];
B1_3_opto = [find(opto1_3.LogicalVecs{opton},1), find(opto1_3.LogicalVecs{opton},1,'last')];
B2_3_opto = [find(opto2_3.LogicalVecs{opton},1), find(opto2_3.LogicalVecs{opton},1,'last')];


mouse1_t = [sum(B1_1_log(end-5399:end))/5400, sum(B1_2_log(B1_2_opto(1):B1_2_opto(2)))/(B1_2_opto(2)-B1_2_opto(1)+1), sum(B1_2_log(end-5399:end))/5400, sum(B1_3_log(B1_3_opto(1):B1_3_opto(2)))/(B1_3_opto(2)-B1_3_opto(1)+1)];
mouse2_t = [sum(B2_1_log(end-5399:end))/5400, sum(B2_2_log(B2_2_opto(1):B2_2_opto(2)))/(B2_2_opto(2)-B2_2_opto(1)+1), sum(B2_2_log(end-5399:end))/5400, sum(B2_3_log(B2_3_opto(1):B2_3_opto(2)))/(B2_3_opto(2)-B2_3_opto(1)+1)];
