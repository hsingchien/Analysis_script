%% new --> old

sfpid = [];

for i = 1:size(new_centroid,1)
   this_c = new_centroid(i,:);
   dist_to_old = old_centroid - this_c;
   dist_to_old = sqrt(dist_to_old(:,1).^2 + dist_to_old(:,2).^2);
   [~, ids] = sort(dist_to_old);
   sfpid = [sfpid; ids(1)];
   
   
end
ms_new.cell_label= ms_old.cell_label(sfpid);

%% old --> new
old_goodid = find(ms_old.cell_label ==1);
new_centroid = ms_new.centroids_xz;
ms_new.cell_label = ms_new.cell_label * 0;
for i = 1:length(old_goodid)
   old_coord = ms_old.centroids_xz(old_goodid(i),:);
   dist_to_new = new_centroid - old_coord;
   dist_to_new = sqrt(dist_to_new(:,1).^2 + dist_to_new(:,2).^2);
   [~,ids] = sort(dist_to_new);
   ms_new.cell_label(ids(1)) = 1;
end
%%
ms = ms_new;
save('ms_cleaned.mat','ms');