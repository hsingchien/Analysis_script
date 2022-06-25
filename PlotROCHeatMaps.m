% get list of mat files
flist = {dir('*.mat').name};
targetFolder = 'Tuning\';
pthreshold = 97.5;
all_active = []; all_inhib = [];

GenTracker = {}; ATracker = {};
%% load and plot
for i = 1:numel(flist)
    load(flist{i});
    M1 = rocStruct{1}; M2 = rocStruct{2};
    if isempty(M1.ROC)
        continue
    end
    
    if isfield(M1,'GenType')
        m1G = M1.GenType; m2G = M2.GenType;
    else
        m1G = ''; m2G = '';
    end
    session = 'toy';
    self1 = getSigCells(M1.ROC.(session).other.self.Encoding,pthreshold);
    heatMap(self1,'xTickLabels',M1.ROC.(session).other.self.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
    set(gcf,'Position',[0,0,800,700]); saveas(gcf,[targetFolder,'self\',M1.ExperimentID,'-',M1.AnimalID,'-',m1G,'.png'],'png');
    
    self2 = getSigCells(M2.ROC.(session).other.self.Encoding,pthreshold);
    heatMap(self2,'xTickLabels',M2.ROC.(session).other.self.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
    set(gcf,'Position',[0,0,800,700]); saveas(gcf,[targetFolder,'self\',M2.ExperimentID,'-',M2.AnimalID,'-',m2G,'.png'],'png');

    partner1 = getSigCells(M1.ROC.(session).other.partner.Encoding,pthreshold);
    heatMap(partner1,'xTickLabels',M1.ROC.(session).other.partner.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
    set(gcf,'Position',[0,0,800,700]); saveas(gcf,[targetFolder,'partner\',M1.ExperimentID,'-',M1.AnimalID,'-',m1G,'.png'],'png');

    partner2 = getSigCells(M2.ROC.(session).other.partner.Encoding,pthreshold);
    heatMap(partner2,'xTickLabels',M2.ROC.(session).other.partner.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
    set(gcf,'Position',[0,0,800,700]); saveas(gcf,[targetFolder,'partner\',M2.ExperimentID,'-',M2.AnimalID,'-',m2G,'.png'],'png');

    selfdrvd1 = getSigCells(M1.ROC.(session).other.derived_self.Encoding,pthreshold);
    heatMap(selfdrvd1,'xTickLabels',M1.ROC.(session).other.derived_self.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
    set(gcf,'Position',[0,0,600,700]); saveas(gcf,[targetFolder,'self_derived\',M1.ExperimentID,'-',M1.AnimalID,'-',m1G,'.png'],'png');

    selfdrvd2 = getSigCells(M2.ROC.(session).other.derived_self.Encoding,pthreshold);
    heatMap(selfdrvd2,'xTickLabels',M2.ROC.(session).other.derived_self.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
    set(gcf,'Position',[0,0,600,700]); saveas(gcf,[targetFolder,'self_derived\',M2.ExperimentID,'-',M2.AnimalID,'-',m2G,'.png'],'png');

    partnerdrvd1 = getSigCells(M1.ROC.(session).other.derived_partner.Encoding,pthreshold);
    heatMap(partnerdrvd1,'xTickLabels',M1.ROC.(session).other.derived_partner.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
    set(gcf,'Position',[0,0,600,700]); saveas(gcf,[targetFolder,'partner_derived\',M1.ExperimentID,'-',M1.AnimalID,'-',m1G,'.png'],'png');
    
    partnerdrvd2 = getSigCells(M2.ROC.(session).other.derived_partner.Encoding,pthreshold);
    heatMap(partnerdrvd2,'xTickLabels',M2.ROC.(session).other.derived_partner.Behavior.EventNames,'colormap',[1,0,0;1,1,1;0,1,1],'fontSize',10);
    set(gcf,'Position',[0,0,600,700]); saveas(gcf,[targetFolder,'partner_derived\',M2.ExperimentID,'-',M2.AnimalID,'-',m2G,'.png'],'png');
    
    close all
    M1_active = sum(selfdrvd1==1,1)/size(selfdrvd1,1);
    M1_inhib = sum(selfdrvd1==-1,1)/size(selfdrvd1,1);
    M2_active = sum(selfdrvd2==1,1)/size(selfdrvd2,1);
    M2_inhib = sum(selfdrvd2==-1,1)/size(selfdrvd2,1);
    all_active = [all_active; M1_active; M2_active];
    all_inhib = [all_inhib; M1_inhib; M2_inhib];
    GenTracker = [GenTracker, m1G, m2G];
    ExID = strsplit(M1.ExperimentID,'_');
    ATracker = [ATracker, [ExID{1},'-',M1.AnimalID], [ExID{1},'-',M2.AnimalID]];
    

end

