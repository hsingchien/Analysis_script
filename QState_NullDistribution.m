% superimpose individual session

%% circshift null distribution
p_diag = []; p_same = []; p_neighbor = [];
for i = 1:1000
    XZ98qseries_shifted = circshift(XZ98qseries, randi(length(XZ98qseries)));
    q_diff = diff([XZ93qseries,XZ98qseries_shifted],[],2);
    p_time_spent_together = sum(q_diff == 0)/length(q_diff)
    p_time_spent_neighbor = sum(or(abs(q_diff)==1,abs(q_diff)==3))/length(q_diff)
    p_time_spent_diag = sum(abs(q_diff)==2)/length(q_diff)
    p_diag = [p_diag; p_time_spent_diag];
    p_same = [p_same; p_time_spent_together];
    p_neighbor = [p_neighbor; p_time_spent_neighbor];
end
f=figure;ax=axes('Parent',f, 'NextPlot','add');
colors = [0.5,0.5,0;0,0.7,0.4;0,0,0.8];
XZBoxPlot({p_same,p_neighbor,p_diag},[1,2,3],[],{'Same','Neighbor','Diag'},colors,ax);
p1 = plot(ax,[1,2,3], [0.15,0.52,0.32],'Marker','o','MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none','MarkerSize',8,'LineStyle','none');
ylabel('Fraction of time');
%%
p_diag = []; p_same = []; p_neighbor = [];
for i = 1:1000
    XZ98qseries_shifted = circshift(pairqseries(:,1), randi(length(pairqseries)));
    q_diff = diff([XZ98qseries_shifted,pairqseries(:,2)],[],2);
    p_time_spent_together = sum(q_diff == 0)/length(q_diff);
    p_time_spent_neighbor = sum(or(abs(q_diff)==1,abs(q_diff)==3))/length(q_diff);
    p_time_spent_diag = sum(abs(q_diff)==2)/length(q_diff);
    p_diag = [p_diag; p_time_spent_diag];
    p_same = [p_same; p_time_spent_together];
    p_neighbor = [p_neighbor; p_time_spent_neighbor];
end
f=figure;ax=axes('Parent',f, 'NextPlot','add');
colors = [0.5,0.5,0;0,0.7,0.4;0,0,0.8];
XZBoxPlot({p_same,p_neighbor,p_diag},[1,2,3],[],{'Same','Neighbor','Diag'},colors,ax);
p1 = plot(ax,[1,2,3], [0.15,0.52,0.32],'Marker','o','MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none','MarkerSize',8,'LineStyle','none');
ylabel('Fraction of time');