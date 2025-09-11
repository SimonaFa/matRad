matRad_rc;
matRad_cfg = MatRad_Config.instance();
% matRad_cfg.propOpt.defaultMaxIter = 500000;
% matRad_cfg.propOpt.defaultAccChangeTol = 1e-06;
matRad_cfg.defaults.propOpt.maxIter = 500000;

load('PROSTATE.mat');

pln.radiationMode = 'carbon';
pln.machine       = 'Generic_clusterDose_prestep_RBEinterp';
pln.propDoseCalc.calcLET = 0;

pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);

pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 8;
pln.propDoseCalc.doseGrid.resolution.y = 8;
pln.propDoseCalc.doseGrid.resolution.z = 8;

pln.multScen = matRad_multScen(ct, 'nomScen');

%% stf
stf = matRad_generateStf(ct,cst,pln);

%% Set cst
% cst{1,6} = {};
% cst{2,6} = {};
% 
% cst{1,6} = {struct(DoseObjectives.matRad_SquaredOverdosing(5,4.2))};
%cst{2,6}{1}.robustness = 'none';
%cst{3,6}{1}.robustness = 'none';
cst{1,6}{1}.robustness = 'none';

cst{6,6}{1}.robustness = 'none';
cst{7,6}{1}.robustness = 'none';
cst{8,6}{1}.robustness = 'none';
cst{9,6}{1}.robustness = 'none';


cst{1,6}{1}.quantity   = 'physicalDose';
%cst{2,6}{1}.quantity   = 'physicalDose';
%cst{3,6}{1}.quantity   = 'physicalDose';

cst{6,6}{1}.quantity   = 'physicalDose';
cst{7,6}{1}.quantity   = 'physicalDose';
cst{8,6}{1}.quantity   = 'physicalDose';
cst{9,6}{1}.quantity   = 'physicalDose';

% cst{2,6} = {struct(DoseObjectives.matRad_SquaredDeviation(100,42))};
% cst{2,6}{1}.robustness = 'none';
% cst{2,6}{1}.quantity   = 'physicalDose';

%% Dose calc

% Select biological model
% % Available models are:
% %   RBEminMax (LET based): MCN, WED, CAR, LSM, (protons)
% %                          HEL                 (helium)
% %   kernel based:          LEM                 (carbon)
% %   Tabulated              TAB                  (proton, helium, carbon) (requires RBEtable and spectra in base data)
% 
% % Example 1: load the MCN model. The machine input is optional in this case.
% % The output is an instance of the model. 
% pln.bioModel = matRad_bioModel(pln.radiationMode,'MCN', 'generic_MCsquare');
% 
% % Alternatively, the biological model can also be assigned as a structure.
% % The only andatory field in this case is the 'model' field:
% % pln.bioModel.model = 'MCN';
% 
% % Example 2: Some models have parameters that can be tuned by the user.
% % For example, we can instantiate a constRBE model
% %
% % pln.bioModel = matRad_bioModel(pln.radiationMode, 'constRBE');
% %
% % and assign a custom value for the constant RBE
% %
% % pln.bioModel.RBE = 4.2;


% Example 3: Tabulated models require  the additional specification of an
% RBE table and additional base data entries. Additional details are
% provided in the Wiki

pln.bioModel = matRad_bioModel(pln.radiationMode, 'StPo');
%pln.bioModel.fragmentsToInclude  = {'H1','C'};
%pln.bioModel.tableFragmentIndexes  = [1,2,3];
pln.bioModel.RBEtableName       = 'RBEtable_rapidMKM_Kase2008_allIons_MKM_0102';
pln.bioModel.includedFragments = struct('Z',[1,2,3,4,5,6], ...
                                        'A',[1,4,7,9,11,12]);

% ionName = fields(RBEtable.data.alpha);
% RBEtable.data.alpha1 = [];
% for i = 1:size(ionName,1)
%     RBEtable.data.alpha1 = [RBEtable.data.alpha1 RBEtable.data.alpha.(ionName{i})];
% end
 
pln.propOpt.quantityOpt = 'RBExDose';
pln.propDoseCalc.engine = 'HongPB';
%pln.propDoseCalc.numHistoriesPerBeamlet = 1e6;

dij = matRad_calcDoseInfluence(ct,cst,stf,pln);
%% ResultGUI Calculation if there is already a resultGUI
resultGUI_MKM =  matRad_calcDoseForward(ct,cst,stf,pln,resultGUI_newLEMI.w);

%% Fluence opt

resultGUI_newLEMI = matRad_fluenceOptimization(dij,cst,pln,resultGUI_LEMI.w);


%matRadGUI;

%% asd
figure
matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI.effect,3,80);
title("LEMI")
figure
matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI_LEMII.effect,3,80);
title("LEMII")
figure
matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI_LEMIII.effect,3,80);
title("LEMIII")
figure
matRad_plotSliceWrapper(gca(),ct,cst,1,resultGUI_MKM.effect,3,80);
title("MKM")

%% asd
plane = 3;
slice = 40;
cube = resultGUI_newLEMI.RBExDose;
cube2 = resultGUI_fullTOPAS_LEMI.RBExDose;
cube4 = resultGUI_fullTOPAS_LEMII.RBExDose;

%allcubes = [max(resultGUI_fullTopasLEMI.physicalDose(:)) max(resultGUI_fullTopasLEMII.physicalDose(:)) max(resultGUI_fullTopasLEMIII.physicalDose(:)) max(resultGUI_fullTopasMKM.physicalDose(:))];
%allcubesPD = [max(resultGUI.RBExDose(:)) max(resultGUI_newTopasLEMI.RBExDose(:))];% max(resultGUI_LEMIII.effect(:)) max(resultGUI_MKM.effect(:))];
%allcubesMin = [min(resultGUI_fullTopasLEMI.physicalDose(:)) min(resultGUI_fullTopasLEMII.physicalDose(:)) min(resultGUI_fullTopasLEMIII.physicalDose(:)) min(resultGUI_fullTopasMKM.physicalDose(:))];
%allcubesMinPD = [min(resultGUI.RBExDose(:)) min(resultGUI_newTopasLEMI.RBExDose(:))];% min(resultGUI_LEMIII.effect(:)) min(resultGUI_MKM.effect(:))];
%doseWindow = [min(allcubesMin(:)) max(allcubes(:))];
doseWindow = [min(cube4(:)) max(cube4(:))];
figure
subplot(4,2,1)
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('LEMI RBExDose PB')
zoom(2)

doseWindow = [min(cube4(:)) max(cube4(:))];
%doseWindow = [min(cube2(:)) max(cube2(:))];
subplot(4,2,2)
matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('LEMI RBExDose TP')
zoom(2)

cube3 = resultGUI_LEMII.RBExDose;
%doseWindow = [min(allcubesMin(:)) max(allcubes(:))];
doseWindow = [min(cube4(:)) max(cube4(:))];
subplot(4,2,3)
matRad_plotSliceWrapper(gca,ct,cst,1,cube3,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('LEMII RBExDose PB')
zoom(2)

cube4 = resultGUI_fullTOPAS_LEMII.RBExDose;
%doseWindow = [min(allcubesMin(:)) max(allcubes(:))];
doseWindow = [min(cube4(:)) max(cube4(:))];
subplot(4,2,4)
matRad_plotSliceWrapper(gca,ct,cst,1,cube4,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('LEMII RBExDose TP')
zoom(2)

cube5 = resultGUI_LEMIII.RBExDose;
%doseWindow = [min(allcubesMin(:)) max(allcubes(:))];
doseWindow = [min(cube4(:)) max(cube4(:))];
subplot(4,2,5)
matRad_plotSliceWrapper(gca,ct,cst,1,cube5,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('LEMIII RBExDose PB')
zoom(2)

cube6 = resultGUI_fullTOPAS_LEMIII.RBExDose;
%doseWindow = [min(allcubesMin(:)) max(allcubes(:))];
doseWindow = [min(cube4(:)) max(cube4(:))];
subplot(4,2,6)
matRad_plotSliceWrapper(gca,ct,cst,1,cube6,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('LEMIII RBExDose TP')
zoom(2)

cube7 = resultGUI_MKM.RBExDose;
%doseWindow = [min(allcubesMin(:)) max(allcubes(:))];
doseWindow = [min(cube4(:)) max(cube4(:))];
subplot(4,2,7)
matRad_plotSliceWrapper(gca,ct,cst,1,cube7,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('MKM RBExDose PB')
zoom(2)

cube8 = resultGUI_fullTOPAS_MKM.RBExDose;
%doseWindow = [min(allcubesMin(:)) max(allcubes(:))];
doseWindow = [min(cube4(:)) max(cube4(:))];
subplot(4,2,8)
matRad_plotSliceWrapper(gca,ct,cst,1,cube8,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('MKM RBExDose TP')
zoom(2)

%% TOPAS
pln_TOPAS.propDoseCalc.engine = 'TOPAS';

if pln_TOPAS.propDoseCalc.engine == 'TOPAS'
    if pln.bioModel.model == 'StPo'
        pln_TOPAS = pln;
        pln_TOPAS.propDoseCalc.engine = 'TOPAS';
        pln_TOPAS.bioModel = matRad_bioModel(pln.radiationMode, 'TAB');
        pln_TOPAS.bioModel.RBEtableName       = 'RBEtable_rapidMKM_Kase2008_allIons_MKM_0102';
        pln_TOPAS.bioModel.includedFragments = struct('Z',[1,2,3,4,5,6], ...
                                        'A',[1,4,7,9,11,12]);
        pln_TOPAS.bioModel = pln_TOPAS.bioModel.setFragmentIndexesInBaseData(pln.bioModel.fragmentIndexesInBaseData);
        pln_TOPAS.bioModel = pln_TOPAS.bioModel.setTissueAlphaX(pln.bioModel.tissueAlphaX);
        pln_TOPAS.bioModel = pln_TOPAS.bioModel.setTissueBetaX(pln.bioModel.tissueBetaX);

    end
end


pln_TOPAS.propDoseCalc.externalCalculation = 'write';
pln_TOPAS.propDoseCalc.calcBioDose = true;
pln_TOPAS.propDoseCalc.numHistoriesDirect = 100000;

resultGUI_TOPAS_LEMIII = matRad_calcDoseForward(ct,cst,stf,pln_TOPAS_LEMIII,resultGUI_newLEMI.w);

%% if TOPAS calculation is done
pln_TOPAS_MKM.propDoseCalc.externalCalculation = 'C:\Users\l813r\Documents\GitHub\Lisa\userdata\Results\carbon_Generic_clusterDose_prestep_RBEinterp_MKM_Prostate_23-07-25';
resultGUI_fullTOPAS_MKM = matRad_calcDoseForward(ct,cst,stf,pln_TOPAS_MKM,resultGUI_newLEMI.w);

%% 24.06. TG119_LEMI is Box_LEMI --> wrote the wrong phantom
%%
resultGUI_appendMKM = matRad_appendResultGUI(resultGUI_MKM,resultGUI_fullTOPAS_MKM,true,pln_MKM.propDoseCalc.engine);
matRad_compareDose(resultGUI_appendLEMIII.physicalDose, resultGUI_appendLEMIII.(['physicalDose_' pln_LEMIII.propDoseCalc.engine]), ct, cst, [1, 1, 0] , 'off', pln_LEMIII, [2, 2], 3, 'global');

%% DVH
DVH_phys_LEMI = matRad_calcDVH(cst,resultGUI_newLEMI.physicalDose,'cum');
DVH_phys_LEMII = matRad_calcDVH(cst,resultGUI_LEMII.physicalDose,'cum');
DVH_phys_LEMIII = matRad_calcDVH(cst,resultGUI_LEMIII.physicalDose,'cum');
DVH_phys_MKM = matRad_calcDVH(cst,resultGUI_MKM.physicalDose,'cum');

DVH_phys_TOPAS_LEMI = matRad_calcDVH(cst,resultGUI_fullTOPAS_LEMI.physicalDose,'cum');
DVH_phys_TOPAS_LEMII = matRad_calcDVH(cst,resultGUI_fullTOPAS_LEMII.physicalDose,'cum');
DVH_phys_TOPAS_LEMIII = matRad_calcDVH(cst,resultGUI_fullTOPAS_LEMIII.physicalDose,'cum');
DVH_phys_TOPAS_MKM = matRad_calcDVH(cst,resultGUI_fullTOPAS_MKM.physicalDose,'cum');


DVH_RBE_LEMI = matRad_calcDVH(cst,resultGUI_newLEMI.RBExDose,'cum');
DVH_RBE_LEMII = matRad_calcDVH(cst,resultGUI_LEMII.RBExDose,'cum');
DVH_RBE_LEMIII = matRad_calcDVH(cst,resultGUI_LEMIII.RBExDose,'cum');
DVH_RBE_MKM = matRad_calcDVH(cst,resultGUI_MKM.RBExDose,'cum');

DVH_RBE_TOPAS_LEMI = matRad_calcDVH(cst,resultGUI_fullTOPAS_LEMI.RBExDose,'cum');
DVH_RBE_TOPAS_LEMII = matRad_calcDVH(cst,resultGUI_fullTOPAS_LEMII.RBExDose,'cum');
DVH_RBE_TOPAS_LEMIII = matRad_calcDVH(cst,resultGUI_fullTOPAS_LEMIII.RBExDose,'cum');
DVH_RBE_TOPAS_MKM = matRad_calcDVH(cst,resultGUI_fullTOPAS_MKM.RBExDose,'cum');

DVH_BED_LEMI = matRad_calcDVH(cst,resultGUI_newLEMI.BED,'cum');
DVH_BED_LEMII = matRad_calcDVH(cst,resultGUI_LEMII.BED,'cum');
DVH_BED_LEMIII = matRad_calcDVH(cst,resultGUI_LEMIII.BED,'cum');
DVH_BED_MKM = matRad_calcDVH(cst,resultGUI_MKM.BED,'cum');

DVH_BED_TOPAS_LEMI = matRad_calcDVH(cst,resultGUI_fullTOPAS_LEMI.BED,'cum');
DVH_BED_TOPAS_LEMII = matRad_calcDVH(cst,resultGUI_fullTOPAS_LEMII.BED,'cum');
DVH_BED_TOPAS_LEMIII = matRad_calcDVH(cst,resultGUI_fullTOPAS_LEMIII.BED,'cum');
DVH_BED_TOPAS_MKM = matRad_calcDVH(cst,resultGUI_fullTOPAS_MKM.BED,'cum');


%%
figure,
subplot(2,2,1)
vois = [1,6,7,8,9];
for i  = vois
    plot(DVH_phys_LEMI(i).doseGrid,DVH_phys_LEMI(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(DVH_phys_TOPAS_LEMI(i).doseGrid,DVH_phys_TOPAS_LEMI(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
end 
title('LEMI physicalDose')
legend('Pencil Beam','TOPAS')
subplot(2,2,2)
vois = [1,6,7,8,9];
for i  = vois
    plot(DVH_phys_LEMII(i).doseGrid,DVH_phys_LEMII(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(DVH_phys_TOPAS_LEMII(i).doseGrid,DVH_phys_TOPAS_LEMII(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
end 
title('LEMII physicalDose')
legend('Pencil Beam','TOPAS')
subplot(2,2,3)
vois = [1,6,7,8,9];
for i  = vois
    plot(DVH_phys_LEMIII(i).doseGrid,DVH_phys_LEMIII(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(DVH_phys_TOPAS_LEMIII(i).doseGrid,DVH_phys_TOPAS_LEMIII(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
end 
title('LEMIII physicalDose')
legend('Pencil Beam','TOPAS')
subplot(2,2,4)
vois = [1,6,7,8,9];
for i  = vois
    plot(DVH_phys_MKM(i).doseGrid,DVH_phys_MKM(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(DVH_phys_TOPAS_MKM(i).doseGrid,DVH_phys_TOPAS_MKM(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
end 
title('MKM physicalDose')
legend('Pencil Beam','TOPAS')

%% DVH RBExDose

figure,
subplot(2,2,1)
vois = [1,6,7,8,9];
for i  = vois
    plot(DVH_RBE_LEMI(i).doseGrid,DVH_RBE_LEMI(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(DVH_RBE_TOPAS_LEMI(i).doseGrid,DVH_RBE_TOPAS_LEMI(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
end 
hold on 
c =1;
leg = {};
names = {};
%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst{i,2}; 
    c = c+1;

end
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.2);
names{c} = 'Pencil Beam';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.2);
names{c} = 'TOPAS';
c = c+1;
axis([0,12,0,100])
legend ([leg{:}],names ) 
xlabel(' RBExDose [Gy]')
ylabel('Volume [%]')
title('LEMI')
set(gca,'FontSize',14)
grid on


subplot(2,2,2)
vois = [1,6,7,8,9];
for i  = vois
    plot(DVH_RBE_LEMII(i).doseGrid,DVH_RBE_LEMII(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(DVH_RBE_TOPAS_LEMII(i).doseGrid,DVH_RBE_TOPAS_LEMII(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
end 
hold on 
c =1;
leg = {};
names = {};
%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst{i,2}; 
    c = c+1;

end
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.2);
names{c} = 'Pencil Beam';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.2);
names{c} = 'TOPAS';
c = c+1;
axis([0,12,0,100])
legend ([leg{:}],names ) 
xlabel(' RBExDose [Gy]')
ylabel('Volume [%]')
title('LEMII')
set(gca,'FontSize',14)
grid on



subplot(2,2,3)
vois = [1,6,7,8,9];
for i  = vois
    plot(DVH_RBE_LEMIII(i).doseGrid,DVH_RBE_LEMIII(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(DVH_RBE_TOPAS_LEMIII(i).doseGrid,DVH_RBE_TOPAS_LEMIII(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
end 
hold on 
c =1;
leg = {};
names = {};
%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst{i,2}; 
    c = c+1;

end
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.2);
names{c} = 'Pencil Beam';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.2);
names{c} = 'TOPAS';
c = c+1;
axis([0,12,0,100])
legend ([leg{:}],names ) 
xlabel(' RBExDose [Gy]')
ylabel('Volume [%]')
title('LEMIII')
set(gca,'FontSize',14)
grid on


subplot(2,2,4)
vois = [1,6,7,8,9];
for i  = vois
    plot(DVH_RBE_MKM(i).doseGrid,DVH_RBE_MKM(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(DVH_RBE_TOPAS_MKM(i).doseGrid,DVH_RBE_TOPAS_MKM(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
end 
hold on 
c =1;
leg = {};
names = {};
%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst{i,2}; 
    c = c+1;

end
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.2);
names{c} = 'Pencil Beam';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.2);
names{c} = 'TOPAS';
c = c+1;
axis([0,12,0,100])
legend ([leg{:}],names ) 
xlabel(' RBExDose [Gy]')
ylabel('Volume [%]')
title('MKM')
set(gca,'FontSize',14)
grid on

%% DVH BED

figure,
subplot(2,2,1)
vois = [1,6,7,8,9];
for i  = vois
    plot(DVH_BED_LEMI(i).doseGrid,DVH_BED_LEMI(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(DVH_BED_TOPAS_LEMI(i).doseGrid,DVH_BED_TOPAS_LEMI(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
end 
hold on 
c =1;
leg = {};
names = {};
%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst{i,2}; 
    c = c+1;

end
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.2);
names{c} = 'Pencil Beam';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.2);
names{c} = 'TOPAS';
c = c+1;
axis([0,50,0,100])
legend ([leg{:}],names ) 
xlabel(' BED [Gy]')
ylabel('Volume [%]')
title('LEMI')
set(gca,'FontSize',14)
grid on


subplot(2,2,2)
vois = [1,6,7,8,9];
for i  = vois
    plot(DVH_BED_LEMII(i).doseGrid,DVH_BED_LEMII(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(DVH_BED_TOPAS_LEMII(i).doseGrid,DVH_BED_TOPAS_LEMII(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
end 
hold on 
c =1;
leg = {};
names = {};
%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst{i,2}; 
    c = c+1;

end
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.2);
names{c} = 'Pencil Beam';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.2);
names{c} = 'TOPAS';
c = c+1;
axis([0,50,0,100])
legend ([leg{:}],names ) 
xlabel(' BED [Gy]')
ylabel('Volume [%]')
title('LEMII')
set(gca,'FontSize',14)
grid on



subplot(2,2,3)
vois = [1,6,7,8,9];
for i  = vois
    plot(DVH_BED_LEMIII(i).doseGrid,DVH_BED_LEMIII(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(DVH_BED_TOPAS_LEMIII(i).doseGrid,DVH_BED_TOPAS_LEMIII(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
end 
hold on 
c =1;
leg = {};
names = {};
%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst{i,2}; 
    c = c+1;

end
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.2);
names{c} = 'Pencil Beam';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.2);
names{c} = 'TOPAS';
c = c+1;
axis([0,50,0,100])
legend ([leg{:}],names ) 
xlabel(' BED [Gy]')
ylabel('Volume [%]')
title('LEMIII')
set(gca,'FontSize',14)
grid on


subplot(2,2,4)
vois = [1,6,7,8,9];
for i  = vois
    plot(DVH_BED_MKM(i).doseGrid,DVH_BED_MKM(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', '-.');
    hold on
    plot(DVH_BED_TOPAS_MKM(i).doseGrid,DVH_BED_TOPAS_MKM(i).volumePoints,'LineWidth',1.2,'Color',cst{i,5}.visibleColor,'LineStyle', ':');
    hold on
end 
hold on 
c =1;
leg = {};
names = {};
%custom legend
for i = vois
    leg{c} = plot(nan,'Color',cst{i,5}.visibleColor,'LineWidth',1.2);
    hold on 
    names{c} = cst{i,2}; 
    c = c+1;

end
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle','-.','LineWidth',1.2);
names{c} = 'Pencil Beam';
c = c+1;
leg{c} = plot(NaN, 'Color',[0,0,0],'LineStyle',':','LineWidth',1.2);
names{c} = 'TOPAS';
c = c+1;
axis([0,50,0,100])
legend ([leg{:}],names ) 
xlabel(' BED [Gy]')
ylabel('Volume [%]')
title('MKM')
set(gca,'FontSize',14)
grid on

%%
effect =@(n,d,a,b)(n(a*d +b*d*d ));
BED_LEMI = effect(30,resultGUI_newLEMI.RBExDose,0.1,0.05)./0.1 ;