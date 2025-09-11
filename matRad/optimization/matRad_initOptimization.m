function [dij,cst,pln,wInit,optiProb,FLAG_ROB_OPT] = matRad_initOptimization(dij,cst,pln,wInit)
% matRad inverse planning wrapper function
%
% call
%   [dij,cst,pln,wInit,optiProb,FLAG_ROB_OPT] = matRad_initOptimization(dij,cst,pln)
%   [dij,cst,pln,wInit,optiProb,FLAG_ROB_OPT] = matRad_initOptimization(dij,cst,pln,wInit)
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
% Copyright 2024 the matRad development team.
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

% consider VOI priorities
cst  = matRad_setOverlapPriorities(cst);

% check & adjust objectives and constraints internally for fractionation
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
        if ~isa(obj,'matRad_DoseOptimizationFunction')
            try
                obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
            catch
                matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
            end
        end

        if isa(obj, 'matRad_DoseOptimizationFunction')
            obj = obj.setDoseParameters(obj.getDoseParameters()/pln.numOfFractions);
        end
        cst{i,6}{j} = obj;
    end

end

% resizing cst to dose cube resolution
cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,  dij.ctGrid.y,  dij.ctGrid.z,...
    dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

if ~isfield(pln,'bioModel')
    pln.bioModel = 'none';
end

if ~isa(pln.bioModel,'matRad_BiologicalModel')
    pln.bioModel = matRad_BiologicalModel.validate(pln.bioModel,pln.radiationMode);
end


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
            if isa(fObjCell, 'DoseObjectives.matRad_DoseObjectives') || isa(fObjCell, 'DoseConstraints.matRad_DoseConstraint')
                dParams = fObjCell{1}.getDoseParameters();
                %Don't care for Inf constraints
                dParams = dParams(isfinite(dParams));
                %Add do dose list
                fDoses = [fDoses dParams];
            end
        end

        %%%
        if isempty(fDoses)
            fDoses = 2;
        end
        %%%
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

% Instantiate backprojetion
backProjection = matRad_BackProjection();

% Initilaize optimization quantities
[optQuantities, constraintQuantities] = backProjection.getOptimizationConstraintQuantitiesFromCst(cst,pln);
backProjection.instantiateQuatities(optQuantities,constraintQuantities,dij,cst);

% Initilaize weights
if ~exist('wInit', 'var') || isempty(wInit)
    wInit = [];
end

% Check minimum biological quantities available
%for qtIdx=backProjection.quantities'
% For now like this. TODO: loop over the quantities, see what the quantity
% needs. Requires implementation of speific declarations in the quantity
% class. In general could add a validation function in the quantity class
% that checks the dij fields are there. DijDistrib quantitiy already does
% the check in the constructor.
if ~isempty(intersect([backProjection.optimizationQuantities,backProjection.constrainedQuantities], {'effect', 'RBExDose', 'BED'})) && ~all(isfield(dij,{'ax','bx'}))
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

if all(isfield(dij,{'ax','bx'}))
    for i = 1:numel(cst{1,4})
        dij.ixDose{i}  = dij.bx{i}~=0;

        %pre-calculations
        dij.gamma{i}             = zeros(dij.doseGrid.numOfVoxels,dij.numOfScenarios);
        dij.gamma{i}(dij.ixDose{i}) = dij.ax{i}(dij.ixDose{i})./(2*dij.bx{i}(dij.ixDose{i}));

    end
end
% Easier if backprojection already knows about the optimized quantities and
% so on
wInit = backProjection.initializeWeights(cst,dij);
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

optiProb = matRad_OptimizationProblem(backProjection);

if isfield(pln,'propOpt') && isfield(pln.propOpt,'useLogSumExpForRobOpt')
    optiProb.useLogSumExpForRobOpt = pln.propOpt.useLogSumExpForRobOpt;
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
end

end