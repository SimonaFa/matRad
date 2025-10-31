clear;
matRad_rc;
matRad_cfg = MatRad_Config.instance();

load('BOXPHANTOM.mat');

%% Set cst

%cst{1,6} = {struct(DoseObjectives.matRad_SquaredOverdosing(5,4.2))};
%
for i=1:size(cst,1)
    cst{i,6}{1}.robustness = 'none';
    cst{i,6}{1}.quantity   = 'physicalDose';
end
%{
for i=1:size(cst,1)
    cst{i,6}{2}.robustness = 'none';
    cst{i,6}{2}.quantity   = 'clusterDose';
end
%}

%% Set cluster dose cst

for i=1:size(cst,1)
    cst{i,6}{1}.robustness = 'none';
    cst{i,6}{1}.quantity   = 'clusterDose';
end

cst{1,6}{1}.parameters  = {[2.0e16]};
cst{1,6}{1}.penalty     = 1e-27;

cst{2,6}{1}.parameters  = {[3.9e16]};
cst{1,6}{1}.penalty     = 1e-25;

%%

pln.radiationMode   = 'carbon';
pln.machine         = 'Generic_clusterDose_prestep';
pln.multScen        = 'nomScen';
pln.numOfFractions  = 30;

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
pln.propDoseCalc.engine = 'HongPB';

%Optimization Settings
pln.propDoseCalc.calcClusterDose            = 1;
pln.propDoseCalc.calcClusterDoseFromFluence = 0;
pln.propDoseCalc.calcCDScatteringFromDose   = 1;

pln.propDoseCalc.cutOffMethod = 'relative';
pln.propDoseCalc.dosimetricLateralCutOff = 0.98;
pln.propDoseCalc.visBoolLateralCutOff       = 1;

pln.propOpt.quantityOpt = 'physicalDose';%'RBExDose';

%% stf
stf = matRad_generateStf(ct,cst,pln);
stf.machine = pln.machine;

%% Dose calc

% Select biological model
% Available models are:
%   RBEminMax (LET based): MCN, WED, CAR, LSM, (protons)
%                          HEL                 (helium)
%   kernel based:          LEM                 (carbon)

% We will compare the MCN model to constRBE.
% We will plan with the the MCN model.
%pln.bioModel = matRad_MCNamara(); 
% pln.bioModel = matRad_bioModel(pln.radiationMode,'LEM');
%altnerative: pln.bioModel = 'MCN';

pln.bioModel = matRad_bioModel(pln.radiationMode,'doseAveragedTabulatedAlphaBeta');
pln.bioModel.quantityTableName = 'RBE_LEM1_update';
pln.bioModel.includedFragments = struct('Z', {1,2,3,4,5,6}, 'A', {1,4,7,9,11,12});%[1, 1];% {'H1', 'C'};
pln.bioModel.stoppingPowerTableName = 'SPtable'; %'RBEtable_rapidLEMI_testTable.mat';
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Fluence optimization
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Now let's recalculate with the constRBE model
pln.bioModel = matRad_ConstantRBE();
pln.bioModel.RBE = 1.1; %1.1 is standard, this is for illustration

resultGUI_recalc = matRad_calcDoseForward(ct,cst,stf,pln,resultGUI.w);

%% Compare Dose distributions

pln.displayQuantity = 'RBExDose';
matRad_compareDose(resultGUI.RBExDose,resultGUI_recalc.RBExDose,ct,cst, [1, 1, 0] , 'off', pln, [3, 3], 3, 'global');