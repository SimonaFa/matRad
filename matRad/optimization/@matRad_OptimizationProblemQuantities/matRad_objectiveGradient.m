function weightGradient = matRad_objectiveGradient(optiProb,w,dij,cst)
% matRad IPOPT callback: gradient function for inverse planning
% supporting mean dose objectives, EUD objectives, squared overdosage,
% squared underdosage, squared deviation and DVH objectives
%
% call
%   g = matRad_gradFuncWrapper(optiProb,w,dij,cst)
%
% input
%   optiProb: option struct defining the type of optimization
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   cst:     matRad cst struct
%
% output
%   g: gradient of objective function
%
% References
%   [1] http://www.sciencedirect.com/science/article/pii/S0958394701000577
%   [2] http://www.sciencedirect.com/science/article/pii/S0360301601025858
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% get current dose / effect / RBExDose vector

optiProb.BP.compute(dij,w);
d = optiProb.BP.d;
%d = optiProb.BP.GetResult();
gGrad = [];


% get the used scenarios
useScen  = optiProb.BP.scenarios;
scenProb = optiProb.BP.scenarioProb;
useNominalCtScen = optiProb.BP.nominalCtScenarios;

% for optQt=optiProb.BP.optimizationQuantities
%     gGrad.(optQt{1}) = cell(numel(useScen),1);
% end
% retrieve matching 4D scenarios
fullScen      = cell(ndims(d),1);
[fullScen{:}] = ind2sub(size(d),useScen);
contourScen   = fullScen{1};

% if isfield(gGrad, 'physicalDose')
%     gGrad.physicalDose          = cell(size(dij.physicalDose));
%     gGrad.physicalDose(useScen) = {zeros(dij.doseGrid.numOfVoxels,1)};
% end
%For probabilistic optimization
vOmega = 0;

%For COWC
f_COWC = zeros(size(dij.physicalDose));

% compute objective function for every VOI.
for  i = 1:size(cst,1)
   
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
        
        % loop over the number of constraints and objectives for the current VOI
        for j = 1:numel(cst{i,6})
            
            %Get current optimization function
            objective = cst{i,6}{j};
            
            % only perform gradient computations for objectives
            if isa(objective,'DoseObjectives.matRad_DoseObjective')
                
                %This is just for temporary compatibility
                % for optQt=optiProb.BP.optimizationQuantities
                %     if any(strcmp(optQt, {'physicalDose', 'RBExD', 'effect', 'BED', 'physicalDoseExp', 'ApproxEffect','MeanAverageEffect'}))
                %         quantityOptimized = optQt{1};
                %     end
                % end
                quantityOptimized = objective.quantity;

                % retrieve the robustness type
                robustness = objective.robustness;
                
                quantityNames = cellfun(@(x) x.quantityName,optiProb.BP.quantities, 'UniformOutput',false);
                quantityOptimizedInstance = optiProb.BP.quantities{strcmp(quantityOptimized,quantityNames)};
                % rescale dose parameters to biological optimization quantity if required
                objective = quantityOptimizedInstance.setBiologicalDosePrescriptions(objective,cst{i,5}.alphaX,cst{i,5}.betaX);
                
                if ~exist('gGrad', 'var') || ~isfield(gGrad,quantityOptimized)
                    if isa(quantityOptimizedInstance, 'matRad_DistributionQuantity')
                        gGrad.(quantityOptimized)          = cell(size(d.(quantityOptimized)));
                        gGrad.(quantityOptimized)(useScen) = {zeros(dij.doseGrid.numOfVoxels,1)};
                    elseif isa(quantityOptimizedInstance, 'matRad_ScalarQuantity')
                        gGrad.(quantityOptimized)                                       = cell(size(d.(quantityOptimized)));
                        gGrad.(quantityOptimized)(:) = {0};           
                    end
                end


                switch robustness
                    case 'none' % if conventional opt: just sum objectiveectives of nominal dose 
                        if isa(quantityOptimizedInstance, 'matRad_DistributionQuantity')
                            for s = useNominalCtScen
                                ixScen = useScen(s);
                                ixContour = contourScen(s);
                                d_i = d.(quantityOptimized){ixScen}(cst{i,4}{ixContour});
                                gGrad.(quantityOptimized){ixScen}(cst{i,4}{ixContour}) = gGrad.(quantityOptimized){ixScen}(cst{i,4}{ixContour}) + objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                            end
                        elseif isa(quantityOptimizedInstance, 'matRad_ScalarQuantity')
                            % Add later the phases
                            d_i = d.(quantityOptimized){i};
                            gGrad.(quantityOptimized){i} = gGrad.(quantityOptimized){i} + objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                         end
                            %add to dose gradient
                            %mean(objective.penalty*objective.computeDoseObjectiveGradient(d_i))
                            %doseGradient{ixScen}(cst{i,4}{ixContour}) = gGrad.physicalDose{ixScen}(cst{i,4}{ixContour}) + objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                    case 'STOCH' % perform stochastic optimization with weighted / random scenarios
                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            
                            d_i = d.(quantityOptimized){ixScen}(cst{i,4}{ixContour});
                            
                            gGrad.(quantityOptimized){ixScen}(cst{i,4}{ixContour}) = gGrad.(quantityOptimized){ixScen}(cst{i,4}{ixContour}) + ...
                                (objective.penalty*objective.computeDoseObjectiveGradient(d_i) * scenProb(s));
                            
                        end
                        
                    case 'PROB' % use the expectation value and the integral variance influence matrix
                        %First check the speficic cache for probabilistic
                        % if ~exist('doseGradientExp','var')
                        %     doseGradientExp{1} = zeros(dij.doseGrid.numOfVoxels,1);
                        % end
                        nPhases = size(d.(quantityOptimized),2);

                        if isa(quantityOptimizedInstance, 'matRad_DistributionQuantity')
                            
                            if nPhases==1
                                structIdxs = cat(1,cst{i,4}{useNominalCtScen});
                                structIdxs = {unique(structIdxs)};
                            else
                                structIdxs = cst{i,4}(useNominalCtScen);
                            end
                        
                            for phaseIdx=1:nPhases

                                d_i = d.(quantityOptimized){1,phaseIdx}(structIdxs{phaseIdx});
        
                                gGrad.(quantityOptimized){phaseIdx}(cst{i,4}{phaseIdx}) = gGrad.(quantityOptimized){phaseIdx}(cst{i,4}{phaseIdx}) + objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                            end
                        
                        else
                            for phaseIdx=1:nPhases
                               d_i = d.(quantityOptimized){i};
                               gGrad.(quantityOptimized){i,phaseIdx} = gGrad.(quantityOptimized){i, phaseIdx} + objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                            end
                        end
                        
                        % 
                        % for phaseIdx=1:nPhases
                        % 
                        %     d_i = d.(quantityOptimized){1,phaseIdx}(structIdxs{phaseIdx});
                        % 
                        %     gGrad.(quantityOptimized){phaseIdx}(cst{i,4}{phaseIdx}) = gGrad.(quantityOptimized){phaseIdx}(cst{i,4}{phaseIdx}) + objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                        % end

                    case 'VWWC'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                        contourIx = unique(contourScen);
                        if ~isscalar(contourIx)
                            % voxels need to be tracked through the 4D CT,
                            % not yet implemented
                            matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                        end
                        
                        % prepare min/max dose vector for voxel-wise worst case
                        if ~exist('d_tmp','var')
                            d_tmp = [d{useScen}];
                        end
                        
                        d_Scen = d_tmp(cst{i,4}{contourIx},:);
                        [d_max,max_ix] = max(d_Scen,[],2);
                        [d_min,min_ix] = min(d_Scen,[],2);
                        
                        if isequal(cst{i,3},'OAR')
                            d_i = d_max;
                        elseif isequal(cst{i,3},'TARGET')
                            d_i = d_min;
                        end
                        
                        if any(isnan(d_i))
                            matRad_cfg.dispWarning('%d NaN values in gradient.',numel(isnan(d_i)));
                        end
                        
                        deltaTmp = objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                        
                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            
                            if isequal(cst{i,3},'OAR')
                                currWcIx = double(max_ix == s);
                            elseif isequal(cst{i,3},'TARGET')
                                currWcIx = double(min_ix == s);
                            end
                            
                            doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + deltaTmp.*currWcIx;
                        end
                        
                    case 'VWWC_INV'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                        contourIx = unique(contourScen);
                        if ~isscalar(contourIx)
                            % voxels need to be tracked through the 4D CT,
                            % not yet implemented
                            matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                        end
                        
                        % prepare min/max dose vector for voxel-wise worst case
                        if ~exist('d_tmp','var')
                            d_tmp = [d{useScen}];
                        end
                        
                        d_Scen = d_tmp(cst{i,4}{1},:);
                        [d_max,max_ix] = max(d_Scen,[],2);
                        [d_min,min_ix] = min(d_Scen,[],2);
                        
                        if isequal(cst{i,3},'OAR')
                            d_i = d_min;
                        elseif isequal(cst{i,3},'TARGET')
                            d_i = d_max;
                        end
                        
                        if any(isnan(d_i))
                            matRad_cfg.dispWarning('%d NaN values in gradFuncWrapper.',numel(isnan(d_i)));
                        end
                        
                        deltaTmp = objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                        
                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            
                            if isequal(cst{i,3},'OAR')
                                currWcIx = double(min_ix == s);
                            elseif isequal(cst{i,3},'TARGET')
                                currWcIx = double(max_ix == s);
                            end
                            
                            doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + deltaTmp.*currWcIx;
                        end
                        
                    case 'COWC' % composite worst case consideres ovarall the worst objective function value
                        %First check the speficic cache for COWC
                        if ~exist('delta_COWC','var')
                            delta_COWC         = cell(size(doseGradient));
                            delta_COWC(useScen)    = {zeros(dij.doseGrid.numOfVoxels,1)};
                        end
                        
                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            
                            d_i = d{ixScen}(cst{i,4}{ixContour});
                            
                            f_COWC(ixScen) = f_COWC(ixScen) + objective.penalty*objective.computeDoseObjectiveFunction(d_i);
                            delta_COWC{ixScen}(cst{i,4}{ixContour}) = delta_COWC{ixScen}(cst{i,4}{ixContour}) + objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                        end
                        
                    case 'OWC' % objective-wise worst case consideres the worst individual objective function value
                        %First check the speficic cache for COWC
                        f_OWC = zeros(size(doseGradient));
                        
                        if ~exist('delta_OWC','var')
                            delta_OWC = cell(size(doseGradient));
                            delta_OWC(useScen) = {zeros(dij.doseGrid.numOfVoxels,1)};
                        end
                        
                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            
                            d_i = d{ixScen}(cst{i,4}{ixContour});
                            
                            f_OWC(ixScen) = objective.penalty*objective.computeDoseObjectiveFunction(d_i);
                            
                            delta_OWC{ixScen}(cst{i,4}{ixContour}) = objective.penalty*objective.computeDoseObjectiveGradient(d_i);
                            
                        end
                          
                        switch optiProb.useMaxApprox
                            case 'logsumexp'
                                [~,fGrad] = optiProb.logSumExp(f_OWC);
                            case 'pnorm'
                                [~,fGrad] = optiProb.pNorm(f_OWC,numel(useScen));
                            case 'none'
                                [~,ix] = max(f_OWC(:));
                                fGrad = zeros(size(f_OWC));
                                fGrad(ix) = 1;
                            case 'otherwise'
                                matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
                                [~,ix] = max(f_OWC(:));
                                fGrad = zeros(size(f_OWC));
                                fGrad(ix) = 1;
                        end
                        
                        for s = 1:numel(useScen)
                            ixScen = useScen(s);
                            ixContour = contourScen(s);
                            if fGrad(ixScen ) ~= 0
                                doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + fGrad(ixScen)*delta_OWC{ixScen}(cst{i,4}{ixContour});
                            end
                        end
                        
                    otherwise
                        matRad_cfg.dispError('Robustness setting %s not supported!',objective.robustness);
                        
                end  %robustness type                              
            elseif isa(objective, 'OmegaObjectives.matRad_OmegaObjective')
            
                robustness = objective.robustness;

                %This is just for temporary compatibility
                % for optQt=optiProb.BP.optimizationQuantities
                %     if any(strcmp(optQt, {'meanVariance', 'vTotApprox', 'vAlpha', 'vBeta', 'MeanAverageEffectVariance', 'AlphaOnlyVariance', 'SqrtBetaOnlyVariance'}))
                %         quantityOptimizedVariance = optQt{1};
                %     end
                % end
                quantityOptimizedVariance = objective.quantity;
                quantityNames = cellfun(@(x) x.quantityName,optiProb.BP.quantities, 'UniformOutput',false);
                quantityOptimizedInstance = optiProb.BP.quantities{strcmp(quantityOptimizedVariance,quantityNames)};
                % rescale dose parameters to biological optimization quantity if required
                %objective = quantityOptimizedInstance.setBiologicalDosePrescriptions(objective,cst{i,5}.alphaX,cst{i,5}.betaX);


                if ~exist('gGrad', 'var') || ~isfield(gGrad,quantityOptimizedVariance)
                    gGrad.(quantityOptimizedVariance)          = cell(size(d.(quantityOptimizedVariance)));
                    gGrad.(quantityOptimizedVariance)(:) = {0};

                end

                 nPhasesOmega = size(d.(quantityOptimizedVariance),2);
                 switch robustness
                     case 'PROB'

                        if nPhasesOmega==1
                            structIdxs = cat(1,cst{i,4}{useNominalCtScen});
                            structIdxs = {unique(structIdxs)};
                        else
                            structIdxs = cst{i,4}(useNominalCtScen);
                        end

                        for phaseIdx = 1:nPhasesOmega
                            vTot = d.(quantityOptimizedVariance){i, phaseIdx};
                            gGrad.(quantityOptimizedVariance){i, phaseIdx} = gGrad.(quantityOptimizedVariance){i,phaseIdx} + objective.penalty * objective.computeTotalVarianceGradient(vTot,structIdxs{phaseIdx});
                        end
                 end

            end
        end %objective loop       
    end %empty check    
end %cst structure loop


if exist('delta_COWC','var')   
    switch optiProb.useMaxApprox
        case 'logsumexp'
            [~,fGrad] = optiProb.logSumExp(f_COWC);
        case 'pnorm'
            [~,fGrad] = optiProb.pNorm(f_COWC,numel(useScen));
        case 'none'
            [~,ixCurrWC] = max(f_COWC(:));
            fGrad = zeros(size(f_COWC));
            fGrad(ixCurrWC) = 1;
        case 'otherwise'
            matRad_cfg.dispWarning('Unknown maximum approximation desired. Using ''none'' instead.');
            [~,ixCurrWC] = max(f_COWC(:));
            fGrad = zeros(size(f_COWC));
            fGrad(ixCurrWC) = 1;
    end
    
    for s = 1:numel(useScen)
        ixScen = useScen(s);
        if fGrad(ixScen) ~= 0
            doseGradient{ixScen} = doseGradient{ixScen} + fGrad(ixScen)*delta_COWC{ixScen};
        end
    end
end

weightGradient = zeros(dij.totalNumOfBixels,1);

% fGrad has an entry for each quantity on which objective functions are
% defined. the fGrad is the sum over all the objfunctions defined for that
% quantity, then gradient wrt the quantity is defined by quantity itself.
% Later on, placing the right fGrad in the right place should be
% determined by the objective itself. Now done manually
%gGrad(find(strcmp(optiProb.BP.optimizationQuantities,quantityOptimized))) = {doseGradient};
optiProb.BP.computeGradient(dij,gGrad,w);
g = optiProb.BP.wGrad;

% for s = 1:numel(useScen)
%     % This could be moved to BP
%     for quantityIdx=optiProb.BP.optimizationQuantities
%         weightGradient = weightGradient + g.(quantityIdx{1}){useScen(s)};
%     end
% end
for qtIdx=optiProb.BP.optimizationQuantities
    nScensOrStructs   = find(cellfun(@(x) ~isempty(x), g.(qtIdx{1})))';
    for elementIdx=nScensOrStructs
        weightGradient = weightGradient + g.(qtIdx{1}){elementIdx};
    end
end


gradientChecker = 0;
if gradientChecker == 1
    f =  matRad_objectiveFunction(optiProb,w,dij,cst);
    epsilon = 1e-3;


    ix = unique(randi([dij.totalNumOfBixels],1,5));

    for i=ix

        wInit = w;
        wInit(i) = wInit(i) + epsilon;
        fDel= matRad_objectiveFunction(optiProb,wInit,dij,cst);
        numGrad = (fDel - f)/epsilon;
        diff = (numGrad/weightGradient(i) - 1)*100;
        fprintf(['grad val #' num2str(i) '- rel diff numerical and analytical gradient = ' num2str(diff) '\n']);
        %fprintf([' any nan or zero for photons' num2str(sum(isnan(glog{1}))) ',' num2str(sum(~logical(glog{1}))) ' for protons: ' num2str(sum(isnan(glog{2}))) ',' num2str(sum(~logical(glog{2}))) '\n']);
    end
end

