function resultGUI = matRad_calcCubes(w,dij,scenNum,boolInterpolate, calcChosenQuantities)
% matRad computation of all cubes for the resultGUI struct
% which is used as result container and for visualization in matRad's GUI
%
% call
%   resultGUI = matRad_calcCubes(w,dij)
%   resultGUI = matRad_calcCubes(w,dij,scenNum)
%
% input
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   scenNum: optional: number of scenario to calculated (default 1)
%
% output
%   resultGUI: matRad result struct
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

if  ~exist('scenNum', 'var') || isempty(scenNum)
    scenNum = 1;
end

if ~exist('boolInterpolate', 'var') || isempty(boolInterpolate)
    boolInterpolate = true;
end

if ~exist('calcChosenQuantities', 'var') || isempty(calcChosenQuantities)
    calcChosenQuantities = matRad_BackProjectionQuantity.getAvailableOptimizationQuantities();
else
    availableQuantities = matRad_BackProjectionQuantity.getAvailableOptimizationQuantities();
    %tmp = struct();
    % Fieldnames
    %fields = fieldnames(availableStruct);
    % Create empty struct with cell2struct
    %tmp = cell2struct(cell(length(fieldnames(availableQuantities)), 1), fieldnames(availableQuantities), 1);
    tmp = [];
    for idx=1:length(calcChosenQuantities)
        if any(strcmp(calcChosenQuantities{idx}, {availableQuantities.quantityName}))
            tmp = [tmp availableQuantities(strcmp(calcChosenQuantities{idx}, {availableQuantities.quantityName}))];
            %tmp(end+1) = availableQuantities(strcmp(calcChosenQuantities{idx}, {availableQuantities.quantityName})); 
        end
    end
    calcChosenQuantities = tmp;
end

resultGUI.w = w;

if isfield(dij,'numParticlesPerMU')
    resultGUI.MU = (w.*1e6) ./ dij.numParticlesPerMU;
end

% get bixel - beam correspondence  
for i = 1:dij.numOfBeams
    beamInfo(i).suffix = ['_beam', num2str(i)];
    beamInfo(i).logIx  = (dij.beamNum == i);
end
beamInfo(dij.numOfBeams+1).suffix = '';
beamInfo(dij.numOfBeams+1).logIx  = true(size(resultGUI.w,1),1);

[ctScen,~] = ind2sub(size(dij.physicalDose),scenNum);

%% Calc single quantities

for idxQ = 1:length(calcChosenQuantities)
    %if strcmp(calcChosenQuantities(idxQ).quantityName, 'physicalDose') ||  strcmp(calcChosenQuantities(idxQ).quantityName, 'constantRBExDose')
    if strcmp(calcChosenQuantities(idxQ).quantityName, 'physicalDose') ||  strcmp(calcChosenQuantities(idxQ).quantityName, 'RBExDose')
        for j=1:length(beamInfo)
            fieldRG = [calcChosenQuantities(idxQ).quantityName, beamInfo(j).suffix];
            cubeQuantity = matRad_calcSingleQuantity(dij,w.*beamInfo(j).logIx,calcChosenQuantities(idxQ),scenNum);
            %resultGUI = matRad_calcSingleQuantity(dij, resultGUI, calcChosenQuantities(idxQ), scenNum, beamInfo);
            resultGUI.([fieldRG]) = cubeQuantity;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Physical Dose
%{
doseFields = {'physicalDose','doseToWater'};
doseQuantities = {'','_std','_batchStd'};
% compute physical dose for all beams individually and together
for j = 1:length(doseFields)
    for k = 1:length(doseQuantities)
        % Check if combination is a field in dij, otherwise skip
        if isfield(dij,[doseFields{j} doseQuantities{k}])
            % Handle standard deviation fields and add quadratically
            if ~isempty(strfind(lower(doseQuantities{1}),'std'))
                for i = 1:length(beamInfo)
                    resultGUI.([doseFields{j}, doseQuantities{k}, beamInfo(i).suffix]) = sqrt(reshape(full(dij.([doseFields{j} doseQuantities{k}]){scenNum}.^2 * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions));
                    resultGUI.([doseFields{j}, doseQuantities{k}, beamInfo(i).suffix])(isnan(resultGUI.([doseFields{j}, doseQuantities{k}, beamInfo(i).suffix]))) = 0;
                end
            % Handle normal fields as usual
            else
                for i = 1:length(beamInfo)
                    resultGUI.([doseFields{j}, doseQuantities{k}, beamInfo(i).suffix]) = reshape(full(dij.([doseFields{j} doseQuantities{k}]){scenNum} * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions);
                end
            end
        end
    end
end
%}

% FIX THIS
%{
if ~isfield(dij,'doseWeightingThreshold')
    dij.doseWeightingThreshold = 0.01;
end
absoluteDoseWeightingThreshold = dij.doseWeightingThreshold*max(resultGUI.physicalDose(:));
%}


%% LET
% consider LET
if isfield(dij,'mLETDose')
    for i = 1:length(beamInfo)
        LETDoseCube                                 = reshape(full(dij.mLETDose{scenNum} * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions);
        resultGUI.(['LET', beamInfo(i).suffix])     = zeros(dij.doseGrid.dimensions);
        ix                                          = resultGUI.(['physicalDose', beamInfo(i).suffix]) > absoluteDoseWeightingThreshold;
        resultGUI.(['LET', beamInfo(i).suffix])(ix) = LETDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix);
    end
end

%% clusterDose
% consider g(I_p)
%{
if isfield(dij,'mClusterDose')
    for i = 1:length(beamInfo)
        clusterDoseCube                                     = reshape(full(dij.mClusterDose{scenNum} * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions);
        resultGUI.(['clusterDose', beamInfo(i).suffix])     = clusterDoseCube;%zeros(dij.doseGrid.dimensions);
        %if isfield(dij, 'doseWeightingThreshold')%.doseWeightingThreshold*max(resultGUI.physicalDose(:));)
        %    ix                                                  = resultGUI.(['clusterDose', beamInfo(i).suffix]) > dij.doseWeightingThreshold*max(resultGUI.clusterDose(:));
        %    resultGUI.(['clusterDose', beamInfo(i).suffix])(ix) = clusterDoseCube(ix)./resultGUI.(['clusterDose', beamInfo(i).suffix])(ix);
        %end
    end
end
%}

%% RBE weighted dose
% consider RBE for protons and skip varRBE calculation
if isfield(dij,'RBE') && isscalar(dij.RBE) && ~ any(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'alpha')), fieldnames(dij)))
    for i = 1:length(beamInfo)
        resultGUI.(['RBExDose', beamInfo(i).suffix]) = resultGUI.(['physicalDose', beamInfo(i).suffix]) * dij.RBE;
    end
elseif any(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'alpha')), fieldnames(dij)))

    if isfield(dij,'RBE') && isscalar(dij.RBE)
        for i = 1:length(beamInfo)
            resultGUI.(['constRBExD', beamInfo(i).suffix]) = resultGUI.(['physicalDose', beamInfo(i).suffix]) * dij.RBE;
        end
    end
    % Load RBE models if MonteCarlo was calculated for multiple models
    if isfield(dij,'RBE_models')
        RBE_model = cell(1,length(dij.RBE_models));
        for i = 1:length(dij.RBE_models)
            RBE_model{i} = ['_' dij.RBE_models{i}];
        end
    else
        RBE_model = {''};
    end

    % Loop through RBE models
    for j = 1:length(RBE_model)
        % Check if combination is a field in dij, otherwise skip
        if isfield(dij,['mAlphaDose' RBE_model{j}])
            for i = 1:length(beamInfo)
                % Get weights of current beam
                wBeam = (resultGUI.w .* beamInfo(i).logIx);

                % consider biological optimization
                ix = dij.bx{ctScen} ~= 0 & resultGUI.(['physicalDose', beamInfo(i).suffix])(:) > 0;
                ixWeighted = dij.bx{ctScen} ~= 0 & resultGUI.(['physicalDose', beamInfo(i).suffix])(:) > absoluteDoseWeightingThreshold;

                % Calculate effect from alpha- and sqrtBetaDose
                resultGUI.(['effect', RBE_model{j}, beamInfo(i).suffix])                = full(dij.(['mAlphaDose' RBE_model{j}]){scenNum} * wBeam + (dij.(['mSqrtBetaDose' RBE_model{j}]){scenNum} * wBeam).^2);
                resultGUI.(['effect', RBE_model{j}, beamInfo(i).suffix])                = reshape(resultGUI.(['effect', RBE_model{j}, beamInfo(i).suffix]),dij.doseGrid.dimensions);

                % Calculate RBExDose from the effect
                resultGUI.(['RBExDose', RBE_model{j}, beamInfo(i).suffix])                 = zeros(size(resultGUI.(['effect', RBE_model{j}, beamInfo(i).suffix])));
                resultGUI.(['RBExDose', RBE_model{j}, beamInfo(i).suffix])(ix)             = (sqrt(dij.ax{ctScen}(ix).^2 + 4 .* dij.bx{ctScen}(ix) .* resultGUI.(['effect', RBE_model{j}, beamInfo(i).suffix])(ix)) - dij.ax{ctScen}(ix))./(2.*dij.bx{ctScen}(ix));

                % Divide RBExDose with the physicalDose to get the plain RBE cube
                resultGUI.(['RBE', RBE_model{j}, beamInfo(i).suffix])                   = zeros(size(resultGUI.(['effect', RBE_model{j}, beamInfo(i).suffix])));
                resultGUI.(['RBE', RBE_model{j}, beamInfo(i).suffix])(ixWeighted)               = resultGUI.(['RBExDose', RBE_model{j}, beamInfo(i).suffix])(ixWeighted)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ixWeighted);

                % Initialize alpha/beta cubes
                resultGUI.(['alpha', RBE_model{j}, beamInfo(i).suffix])                 = zeros(dij.doseGrid.dimensions);
                resultGUI.(['beta',  RBE_model{j}, beamInfo(i).suffix])                 = zeros(dij.doseGrid.dimensions);
                resultGUI.(['alphaDoseCube', RBE_model{j}, beamInfo(i).suffix])         = zeros(dij.doseGrid.dimensions);
                resultGUI.(['SqrtBetaDoseCube',  RBE_model{j}, beamInfo(i).suffix])     = zeros(dij.doseGrid.dimensions);

                % Calculate alpha and weighted alphaDose
                AlphaDoseCube                                                           = full(dij.(['mAlphaDose' RBE_model{j}]){scenNum} * wBeam);
                resultGUI.(['alpha', RBE_model{j}, beamInfo(i).suffix])(ix)             = AlphaDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix);
                resultGUI.(['alphaDoseCube', RBE_model{j}, beamInfo(i).suffix])(ix)     = AlphaDoseCube(ix);

                % Calculate beta and weighted sqrtBetaDose
                SqrtBetaDoseCube                                                        = full(dij.(['mSqrtBetaDose' RBE_model{j}]){scenNum} * wBeam);
                resultGUI.(['beta', RBE_model{j}, beamInfo(i).suffix])(ix)              = (SqrtBetaDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix)).^2;
                resultGUI.(['SqrtBetaDoseCube', RBE_model{j}, beamInfo(i).suffix])(ix)  = SqrtBetaDoseCube(ix);
            end
        end
    end
end

%% Calculate Biological Effective Dose (BED)

% When depth Dependent alpha beta values are calculated in dij calculation
if isfield(dij,'ax') && isfield(dij,'bx')
    ixWeighted = dij.ax{ctScen} > 0 & dij.bx{ctScen} > 0 & resultGUI.(['physicalDose', beamInfo(i).suffix])(:) > absoluteDoseWeightingThreshold;

    if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
        for i = 1:length(beamInfo)
            % photon equivaluent BED = n * effect / alphax
            resultGUI.(['BED', beamInfo(i).suffix]) = zeros(dij.doseGrid.dimensions);
            resultGUI.(['BED', beamInfo(i).suffix])(ixWeighted) = full(resultGUI.(['effect', beamInfo(i).suffix])(ixWeighted) ./dij.ax{ctScen}(ixWeighted));
            resultGUI.(['BED', beamInfo(i).suffix]) = reshape(resultGUI.(['BED', beamInfo(i).suffix]), dij.doseGrid.dimensions);
        end
        matRad_cfg.dispWarning('Photon Equiavlent BED calculated');
    else
        % Get Alpha and Beta Values form dij.ax and dij.bx
        for i = 1:length(beamInfo)
            %         ix = ~isnan(dij.ax{1}./dij.bx{1});
            if isfield(resultGUI, 'RBExDose')
                dose = resultGUI.(['RBExDose', beamInfo(i).suffix]);
            else
                dose = resultGUI.(['physicalDose', beamInfo(i).suffix]);
            end
            effect = dij.ax{ctScen}.* dose(:) + dij.bx{ctScen}.*dose(:).^2;
            resultGUI.(['BED', beamInfo(i).suffix]) = zeros(dij.doseGrid.dimensions);
            resultGUI.(['BED', beamInfo(i).suffix])(ixWeighted) = effect(ixWeighted)./dij.ax{ctScen}(ixWeighted);
            resultGUI.(['BED', beamInfo(i).suffix]) = reshape(resultGUI.(['BED', beamInfo(i).suffix]), dij.doseGrid.dimensions);
        end
        if isfield(resultGUI, 'RBExDose')
            matRad_cfg.dispWarning('Photon Equiavlent BED calculated');
        end
    end
end

%add some dij meta
if isfield(dij,'meta') && isstruct(dij.meta)
    resultGUI.meta = dij.meta;
end


%% Final processing
% Remove suffix for RBExDose if there's only one available
if any(cellfun(@(teststr) ~isempty(strfind(lower(teststr),'alpha')), fieldnames(dij))) && isfield(dij,'RBE_models') && length(dij.RBE_models) == 1
    % Get fieldnames that include the specified RBE model
    fnames = fieldnames(resultGUI);
    fnames = fnames(cellfun(@(teststr) ~isempty(strfind(lower(teststr),lower(dij.RBE_models{1}))), fnames));

    % Rename fields and remove model specifier if there's only one
    for f = 1:length(fnames)
        resultGUI.(erase(fnames{f},['_',dij.RBE_models{1}])) = resultGUI.(fnames{f});
    end

    % Remove old fields
    resultGUI = rmfield(resultGUI,fnames);
end

% group similar fields together
resultGUI = orderfields(resultGUI);

% interpolation if dose grid does not match ct grid
if boolInterpolate
    if isfield(dij,'ctGrid') && any(dij.ctGrid.dimensions~=dij.doseGrid.dimensions)
        myFields = fieldnames(resultGUI);
        for i = 1:numel(myFields)
            if numel(resultGUI.(myFields{i})) == dij.doseGrid.numOfVoxels
    
                % interpolate!
                resultGUI.(myFields{i}) = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
                    resultGUI.(myFields{i}), ...
                    dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);
    
            end
        end
    end
end

end

