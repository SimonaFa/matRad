function [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln,wInit)
% matRad inverse planning wrapper function
%
% call
%   [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln)
%   [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln,wInit)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   pln:        matRad pln struct
%   wInit:      (optional) custom weights to initialize problems
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%   optimizer:  Used Optimizer Object
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

%{
<<<<<<< HEAD
% consider VOI priorities
cst  = matRad_setOverlapPriorities(cst);

% check & adjust objectives and constraints internally for fractionation 
haveDoseOptimizationFunctions = false;
haveClusterDoseOptimizationFunctions = false;

for i = 1:size(cst,1)
    %Compatibility Layer for old objective format
    if isstruct(cst{i,6})
        cst{i,6} = arrayfun(@matRad_DoseOptimizationFunction.convertOldOptimizationStruct,cst{i,6},'UniformOutput',false);
    end
    for j = 1:numel(cst{i,6})

        obj = cst{i,6}{j};

        %In case it is a default saved struct, convert to object
        %Also intrinsically checks that we have a valid optimization
        %objective or constraint function in the end
        if isstruct(obj)
            if strncmp(obj.className,'DoseObjective',13) || strncmp(obj.className,'DoseConstraint',14)
                try
                    obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
                catch
                    matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
                end
            elseif strncmp(obj.className,'ClusterDoseObjective',20)
                try
                    obj = matRad_ClusterDoseOptimizationFunction.createInstanceFromStruct(obj);
                catch
                    matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
                end
                %obj.setClusterDose
            else 
                %error message, or check implementation on branch research/DADR
            end
        end

        if isa(obj,'matRad_DoseOptimizationFunction')
            obj = obj.setDoseParameters(obj.getDoseParameters()/pln.numOfFractions);
            haveDoseOptimizationFunctions = true;
        end
        if isa(obj,'matRad_ClusterDoseOptimizationFunction')
            obj = obj.setDoseParameters(obj.getDoseParameters()/pln.numOfFractions);
            haveClusterDoseOptimizationFunctions = true;
        end

        cst{i,6}{j} = obj;        
    end
end


% resizing cst to dose cube resolution 
cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
                                 dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

% Get rid of voxels that are not interesting for the optimization problem
if ~isfield(pln,'propOpt') || ~isfield(pln.propOpt, 'clearUnusedVoxels')
    pln.propOpt.clearUnusedVoxels = matRad_cfg.defaults.propOpt.clearUnusedVoxels;
end

if pln.propOpt.clearUnusedVoxels
    dij = matRad_clearUnusedVoxelsFromDij(cst, dij);
end



% find target indices and described dose(s) for weight vector
% initialization
V          = [];
doseTarget = [];
ixTarget   = [];

for i = 1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        V = [V;cst{i,4}{1}];

        %Iterate through objectives/constraints
        fDoses = [];
        for fObjCell = cst{i,6}
            dParams = fObjCell{1}.getDoseParameters();
            %Don't care for Inf constraints
            dParams = dParams(isfinite(dParams));
            %Add do dose list
            fDoses = [fDoses dParams];
        end

        doseTarget = [doseTarget fDoses];
        ixTarget   = [ixTarget i*ones(1,length(fDoses))];
    end
end
[doseTarget,i] = max(doseTarget);
ixTarget       = ixTarget(i);
wOnes          = ones(dij.totalNumOfBixels,1);

%Check how to use 4D data
if isfield(pln,'propOpt') && isfield(pln.propOpt,'scen4D')
    scen4D = pln.propOpt.scen4D;
else
    scen4D = 1; %Use only first 4D scenario for optimization
end

% Workaround until future release with consistent data management
totNumCtScen = size(dij.physicalDose,1);

% Validate / Create Scenario model
if ~isfield(pln,'multScen')
    pln.multScen = 'nomScen';
end

if ~isa(pln.multScen,'matRad_ScenarioModel')
    pln.multScen = matRad_ScenarioModel.create(pln.multScen,struct('numOfCtScen',totNumCtScen));
end

if ~isfield(pln,'bioModel')
    pln.bioModel = 'none';
end

if ~isa(pln.bioModel,'matRad_BiologicalModel')
    pln.bioModel = matRad_BiologicalModel.validate(pln.bioModel,pln.radiationMode);
end

%If "all" provided, use all scenarios
if isequal(scen4D,'all')
    scen4D = 1:totNumCtScen;
end

if ~isfield(pln.propOpt, 'quantityOpt') || isempty(pln.propOpt.quantityOpt)
    pln.propOpt.quantityOpt = pln.bioModel.defaultReportQuantity;
    matRad_cfg.dispWarning('quantityOpt was not provided, using quantity suggested by biological model: %s',pln.propOpt.quantityOpt);    
end

% Check optimization quantity
switch pln.propOpt.quantityOpt
    case 'effect'
        backProjection = matRad_EffectProjection;
    case 'RBExDose'
        %Capture special case of constant RBE
        if isa(pln.bioModel,'matRad_ConstantRBE') || (isstruct(pln.bioModel) && strcmp(pln.bioModel.model, 'constRBE'))
            backProjection = matRad_ConstantRBEProjection;
        else
            backProjection = matRad_VariableRBEProjection;
        end
    case 'cluster_Dose'
        backProjection = matRad_ClusterDoseProjection;
    case 'physicalDose'
        backProjection = matRad_DoseProjection;

    case 'BED'
        backProjection = matRad_BEDProjection;
    otherwise
        warning(['Did not recognize biological setting ''' pln.propOpt.quantityOpt '''!\nUsing physical dose optimization!']);
        backProjection = matRad_DoseProjection;
end

% Check minimum biological quantities available
if isa(backProjection,'matRad_EffectProjection') && ~all(isfield(dij,{'ax','bx'}))
    matRad_cfg.dispWarning('Biological optimization requested, but no ax & bx provided in dij. Getting from cst...');

    %First get the voxels where we need it
    validScen = ~cellfun(@isempty,dij.physicalDose);
    d = cellfun(@(D) D*ones(dij.totalNumOfBixels,1),dij.physicalDose(validScen),'UniformOutput',false);
    d = sum(cell2mat(d'),2);
    ixZeroDose = d == 0;

    numOfCtScenarios = numel(cst{1,4});
    for i = 1:numOfCtScenarios
        dij.ax{i} = zeros(dij.doseGrid.numOfVoxels,1);
        dij.bx{i} = zeros(dij.doseGrid.numOfVoxels,1);

        for v = 1:size(cst,1)
            %We already did the overlap stuff so we do not need to care for
            %overlaps here
            dij.ax{i}(cst{v,4}{i}) = cst{v,5}.alphaX;
            dij.bx{i}(cst{v,4}{i}) = cst{v,5}.betaX;
        end

        dij.ax{i}(ixZeroDose) = 0;
        dij.bx{i}(ixZeroDose) = 0;
    end
end


% calculate initial beam intensities wInit
matRad_cfg.dispInfo('Estimating initial weights... ');

=======
>>>>>>> dev_quantities_RBE_tabModels_copyRemo
%}
if exist('wInit','var')
    [dij,cst,pln,wInit,optiProb,FLAG_ROB_OPT] = matRad_initOptimization(dij,cst,pln,wInit);
else

%{
<<<<<<< HEAD
    
    if isfield(dij, 'mClusterDose')
        if ~isempty(dij.mClusterDose)
            bixelWeight =  (doseTarget)/(mean(dij.mClusterDose{1}(V,:)*wOnes));
            wInit       = wOnes * bixelWeight;
        end
    else
        
        bixelWeight =  (doseTarget)/(mean(dij.physicalDose{1}(V,:)*wOnes)); 
        wInit       = wOnes * bixelWeight;
        
    end
    pln.propOpt.bioOptimization = 'none';
    
    matRad_cfg.dispInfo('chosen uniform weight of %f!\n',bixelWeight);
end


%% calculate probabilistic quantities for probabilistic optimization if at least
% one robust objective is defined

linIxDIJ = find(~cellfun(@isempty,dij.physicalDose(scen4D,:,:)))';

%Only select the indexes of the nominal ct Scenarios
linIxDIJ_nominalCT = find(~cellfun(@isempty,dij.physicalDose(scen4D,1,1)))';

FLAG_CALC_PROB = false;
FLAG_ROB_OPT   = false;


for i = 1:size(cst,1)
    for j = 1:numel(cst{i,6})
        if strcmp(cst{i,6}{j}.robustness,'PROB') && numel(linIxDIJ) > 1
            FLAG_CALC_PROB = true;
        end
        if ~strcmp(cst{i,6}{j}.robustness,'none') && numel(linIxDIJ) > 1
            FLAG_ROB_OPT = true;
        end
    end
end

if FLAG_CALC_PROB
    [dij] = matRad_calculateProbabilisticQuantities(dij,cst,pln);
end


% set optimization options
if ~FLAG_ROB_OPT || FLAG_CALC_PROB     % if multiple robust objectives are defined for one structure then remove FLAG_CALC_PROB from the if clause
    ixForOpt = scen4D;
else
    ixForOpt = linIxDIJ;
end

%Give scenarios used for optimization
backProjection.scenarios    = ixForOpt;
backProjection.scenarioProb = pln.multScen.scenProb;
backProjection.nominalCtScenarios = linIxDIJ_nominalCT;
%backProjection.scenDim      = pln.multScen

optiProb = matRad_OptimizationProblem(backProjection);

if isfield(pln,'propOpt') && isfield(pln.propOpt,'useLogSumExpForRobOpt')
    optiProb.useLogSumExpForRobOpt = pln.propOpt.useLogSumExpForRobOpt;
end

if haveClusterDoseOptimizationFunctions && isfield(dij,"mClusterDose")
    
    optiProb.BP_clusterDose = matRad_ClusterDoseProjection;
    
end


%Get Bounds
if ~isfield(pln.propOpt,'boundMU')
    pln.propOpt.boundMU = false;
end


if pln.propOpt.boundMU
    if (isfield(dij,'minMU') || isfield(dij,'maxMU')) && ~isfield(dij,'numParticlesPerMU')
        matRad_cfg.dispWarning('Requested MU bounds but number of particles per MU not set! Bounds will not be enforced and standard [0,Inf] will be used instead!');
    elseif ~isfield(dij,'minMU') && ~isfield(dij,'maxMU')
        matRad_cfg.dispWarning('Requested MU bounds but machine bounds not defined in dij.minMU & dij.maxMU! Bounds will not be enforced and standard [0,Inf] will be used instead!');
    else
        if isfield(dij,'minMU')
            optiProb.minimumW = dij.numParticlesPerMU .* dij.minMU / 1e6;
            matRad_cfg.dispInfo('Using lower MU bounds provided in dij!\n')
        end

        if isfield(dij,'maxMU')
            optiProb.maximumW = dij.numParticlesPerMU .* dij.maxMU / 1e6;
            matRad_cfg.dispInfo('Using upper MU bounds provided in dij!\n')
        end
    end
else
    matRad_cfg.dispInfo('Using standard MU bounds of [0,Inf]!\n')
=======
%}
    [dij,cst,pln,wInit,optiProb,FLAG_ROB_OPT] = matRad_initOptimization(dij,cst,pln);
end

%Dummy
tmpConstRBExDCheck = cellfun(@(quantity) isa(quantity,'matRad_ConstantRBExDose'), optiProb.BP.quantities, 'UniformOutput',false);
if any([tmpConstRBExDCheck{:}]) && ~isfield(dij,'RBE')
    dij.RBE = 1.1;
%>>>>>>> dev_quantities_RBE_tabModels_copyRemo
end

if ~isfield(pln.propOpt,'optimizer')
    %While the default optimizer is IPOPT, we can try to fallback to
    %fmincon in case it does not work for some reason
    if ~matRad_OptimizerIPOPT.IsAvailable()
        pln.propOpt.optimizer = 'fmincon';
    else
        pln.propOpt.optimizer = 'IPOPT';
    end   
end


switch pln.propOpt.optimizer
    case 'IPOPT'
        optimizer = matRad_OptimizerIPOPT;
    case 'fmincon'
        optimizer = matRad_OptimizerFmincon;
    case 'simulannealbnd'
        optimizer = matRad_OptimizerSimulannealbnd;
    otherwise
        warning(['Optimizer ''' pln.propOpt.optimizer ''' not known! Fallback to IPOPT!']);
        optimizer = matRad_OptimizerIPOPT;
end
        
if ~optimizer.IsAvailable()
    matRad_cfg.dispError(['Optimizer ''' pln.propOpt.optimizer ''' not available!']);
end


optimizer = optimizer.optimize(wInit,optiProb,dij,cst);

wOpt = optimizer.wResult;
info = optimizer.resultInfo;

try
    resultGUI = matRad_calcCubes(wOpt,dij);
catch
    matRad_cfg.dispWarning('Unable to compue calcCubes');
end
resultGUI.wUnsequenced = wOpt;
resultGUI.usedOptimizer = optimizer;
resultGUI.info = info;
resultGUI.info.timePerIteration = resultGUI.info.cpu/resultGUI.info.iter;

if ~exist('computeScenarios', 'var') || isempty(computeScenarios)
    computeScenarios = 1;
end

%Robust quantities
try
    if computeScenarios
        if FLAG_ROB_OPT
            if pln.multScen.totNumScen > 1
                for i = 1:pln.multScen.totNumScen
                    scenSubIx = pln.multScen.linearMask(i,:);
                    resultGUItmp = matRad_calcCubes(wOpt,dij,pln.multScen.sub2scenIx(scenSubIx(1),scenSubIx(2),scenSubIx(3)));
                    resultGUI = matRad_appendResultGUI(resultGUI,resultGUItmp,false,sprintf('scen%d',i));
                end
            end
        end
    end
catch
    matRad_cfg.dispWarning('Unable to compute calcCubes');
end
% unblock mex files
clear mex
