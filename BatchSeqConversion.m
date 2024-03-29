behav  = dir('**\BehavCam_0');
behav = unique({behav.folder});
to_remove = [];
for i = 1:length(behav)
if contains(behav{i},'sep')
to_remove = [to_remove, i];
end
if contains(behav{i},'toy')
to_remove = [to_remove, i];
end
end
behav(to_remove) = [];
targdir = 'Y:\Xingjian\Homecage_dual_corrected_videos\Shank3\DLX-DLX\';
%%
behav = {
'E:\MiniscopeData(processed)\NewCage_free_dual\Shank3\DLX-DLX\XZ155_XZ152(m)\2022_04_22\12_59_14_exp\BehavCam_0';
'E:\MiniscopeData(processed)\NewCage_free_dual\Shank3\DLX-DLX\XZ155_XZ154(m)\2022_04_08\17_48_30_exp\BehavCam_0';
'E:\MiniscopeData(processed)\NewCage_free_dual\Shank3\DLX-DLX\XZ155_XZ151(m)\2022_03_31\15_58_11_exp\BehavCam_0';
'E:\MiniscopeData(processed)\NewCage_free_dual\Shank3\DLX-DLX\XZ152_XZ150(m)\2022_03_14\16_20_26_exp\BehavCam_0';
'E:\MiniscopeData(processed)\NewCage_free_dual\Shank3\DLX-DLX\XZ152_XZ146(m)\2022_03_09\15_20_08_exp\BehavCam_0';
}
%%
% expidx=2;
for i = 1:length(behav)
    cd(behav{i});
if exist('0.avi')
    VideoCombine(pwd, 'b',1, true, 'avi');
end
if isempty(dir('*.seq')) & exist('behav_video.avi')
    seqName = 'behavior.seq';
    aviName = 'behav_video.avi';
    seqIo([seqName],'frImgs',struct('codec','png'),'aviName',[aviName]);
end
try
    movefile behav_video.avi behav_video f;
    delete *.avi;
    movefile behav_video behav_video.avi f;
end
% try
%     temp = strsplit(behav{i},'\');
%     
%     mkdir([targdir,'\',temp{5}]);
%     if ~exist([targdir,'\',temp{5},'\',temp{5},'_',temp{6},'_exp.avi'])
%         copyfile('behav_video.avi',[targdir,'\',temp{5},'\',temp{5},'_',temp{6},'_exp.avi']);
%     else
%         copyfile('behav_video.avi',[targdir,'\',temp{5},'\',temp{5},'_',temp{6},'_exp',num2str(expidx),'.avi']);
%         expidx = expidx+1;
%     end
% end
% expidx = 2;
end
%%
behav  = dir('**\');
behav = unique({behav.folder});
for i = 1:length(behav)
    cd(behav{i});

avis = dir("*.avi");
avis = {avis.name};
for j = 1:numel(avis)
    seqName = strrep(avis{j},'avi','seq');
    seqIo([seqName],'frImgs',struct('codec','png'),'aviName',[avis{j}]);
end
end
