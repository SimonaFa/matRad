load('TG119.mat');
matRad_cfg = MatRad_Config.instance();
%% Pln
pln.radiationMode = 'protons';
pln.machine       = 'Generic';

pln.propDoseCalc.calcLET = 0;
pln.propDoseCalc.engine = 'HongPB';

pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [-45,45];
pln.propStf.couchAngles   = [0, 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

quantityOpt   = 'physicalDose';
modelName     = 'none';

pln.bioParam = matRad_bioModel(pln.radiationMode, modelName);

pln.multScen = matRad_multScen(ct,'nomScen');
%pln.multScen.nSamples=2;

pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

%% stf
stf = matRad_generateStf(ct,cst,pln);

%% Set cst
cst{1,6} = {};
cst{3,6} = {};
cst{2,6} = {};

cst{1,6} = {struct(DoseObjectives.matRad_SquaredOverdosing(5,4.2))};
cst{1,6}{1}.robustness = 'none';
cst{1,6}{1}.quantity   = 'physicalDose';

cst{2,6} = {struct(DoseObjectives.matRad_SquaredDeviation(100,42))};
cst{2,6}{1}.robustness = 'none';
cst{2,6}{1}.quantity   = 'physicalDose';

cst{3,6} = {struct(DoseObjectives.matRad_MeanDose(5,0))};
cst{3,6}{1}.robustness = 'none';
cst{3,6}{1}.quantity   = 'physicalDose';

cst{3,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(10,35));
cst{3,6}{2}.robustness = 'none';
cst{3,6}{2}.quantity   = 'physicalDose';
%% Dose calc
pln.propDoseCalc.calcLET = true;
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Optimization

pln.propOpt.quantityOpt = quantityOpt;
resultGUI = matRad_fluenceOptimizationQuantities(dij,cst,pln);

%% Add LET objective

cst{1,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(1,3.2));
cst{1,6}{2}.robustness = 'none';
cst{1,6}{2}.quantity   = 'LETd';

%% Optimization
resultGUI_LET = matRad_fluenceOptimizationQuantities(dij,cst,pln);
%% Compare

figure;
nexttile;
    ax(1) = gca();
    matRad_plotSliceWrapper(ax(1), ct,cst,1,resultGUI.physicalDose,3,80);
    title('Physical Dose','FontSize',20);
    ylabel('pD optimization', 'FontSize',20);

nexttile;
    ax(2) = gca();
    matRad_plotSliceWrapper(ax(2), ct,cst,1,resultGUI.LET,3,80);
    title('LET','FontSize',20);
    ylabel('');

nexttile;
    ax(3) = gca();
    matRad_plotSliceWrapper(ax(3), ct,cst,1,resultGUI_LET.physicalDose,3,80);
    title('');
    ylabel('LET optimization','FontSize',20);

nexttile;
    ax(4) = gca();
    matRad_plotSliceWrapper(ax(4), ct,cst,1,resultGUI_LET.LET,3,80);
    title('');
    ylabel('');

for i=1:4
    set(ax(i), 'XTick', []);
    set(ax(i), 'YTick', []);
    set(ax(i), 'xlabel', []);

end