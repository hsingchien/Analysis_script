f = figure,
a = axes(f);
a.NextPlot = 'add';
xl = line(a, [0,1],[0,0],[0,0],'Color','r');
yl = line(a, [0,0],[0,1],[0,0],'Color','g');
zl = line(a, [0,0],[0,0],[0,1],'Color','b');

vid = VideoReader('behav_video.avi');
vmat = read(vid);

[idx1, idx2, tdiff] = TStampAlign(timeStamps(:,2), headOrientation(:,1));
headOrientation = headOrientation(idx2,:);
vidmat = vmat(:,:,:,idx1);

f2 = figure;
a2 = axes(f2);
imh = imshow(vidmat(:,:,:,1),'Parent',a2);




for i = 1:size(vidmat,4);
    t = headOrientation(i,1);
    q = quaternion(headOrientation(i,2:end));
    x = [1,0,0];
    y = [0,1,0];
    z = [0,0,1];
    xq = rotateframe(q,x);
    yq = rotateframe(q,y);
    zq = rotateframe(q,z);
    set(xl, 'XData', [0, xq(1)], 'YData', [0,xq(2)],'ZData',[0,xq(3)]);
    set(yl, 'XData', [0, yq(1)], 'YData', [0,yq(2)],'ZData',[0,yq(3)]);
    set(zl, 'XData', [0, zq(1)], 'YData', [0,zq(2)],'ZData',[0,zq(3)]);
    set(imh,'CData', vidmat(:,:,:,i));
    drawnow;
    pause(0.03);
end