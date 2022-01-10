behav  = dir('**\BehavCam_0');
behav = unique({behav.folder});
to_remove = [];
for i = 1:length(behav)
if contains(behav{i},'sep')
to_remove = [to_remove, i];
end
end
behav(to_remove) = [];
for i = 1:length(behav)
cd(behav{i});
if ~exist('behav_video.avi')
VideoCombine(pwd, 'b',1, true, 'avi');
end
if isempty(dir('*.seq'))
seqName = 'behavior.seq';
aviName = 'behav_video.avi';
seqIo([seqName],'frImgs',struct('codec','png'),'aviName',[aviName]);
end
movefile behav_video.avi behav_video f;
delete *.avi;
movefile behav_video behav_video.avi f;
end