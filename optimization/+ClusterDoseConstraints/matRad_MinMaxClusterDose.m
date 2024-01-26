classdef matRad_MinMaxClusterDose < ClusterDoseConstraints.matRad_ClusterDoseConstraint
    % matRad_MinMaxDose Implements a MinMaxDose constraint
    %   See matRad_DoseConstraint for interface description
    %
    % use log sum exp approximation, see appendix A in
    % http://scitation.aip.org/content/aapm/journal/medphys/41/8/10.1118/1.4883837
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Constant)
        name = 'Min/Max cluster dose constraint';
        parameterNames = {'cd^{min}', 'cd^{max}','method'};
        parameterTypes = {'clusterdose','clusterdose',{'approx','voxelwise'}};
    end
    
    properties
        parameters = {0,30,1};
        epsilon = 1e-3; %slack parameter for the logistic approximation
    end
    
    methods
        function obj = matRad_MinMaxClusterDose(minClusterDose,maxClusterDose,method)
            
            %If we have a struct in first argument
            if nargin == 1 && isstruct(minClusterDose)
                inputStruct = minDose;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@ClusterDoseConstraints.matRad_ClusterDoseConstraint(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                
                if nargin < 3 || ~ischar(method)
                    method = 'approx';
                end
                
                methodIx = find(strcmp(method,obj.parameterTypes{3}));
                
                if isempty(methodIx) || numel(methodIx) > 1
                    methodIx = 1;
                    msg = ['Dose Constraint method can only be ', strjoin(obj.parameterTypes{3},' or '), '! Using method ''', obj.parameterTypes{3}{methodIx}, '''.'];
                    warning(msg);
                end
                
                obj.parameters{3} = methodIx;
                
                if nargin >= 2 && isscalar(maxClusterDose)
                    obj.parameters{2} = maxClusterDose;
                end
                
                if nargin >= 1 && isscalar(minClusterDose)
                    obj.parameters{1} = minClusterDose;
                end
            end
        end
        
        %Overloads the struct function to add constraint specific
        %parameters
        function s = struct(obj)
            s = struct@ClusterDoseConstraints.matRad_ClusterDoseConstraint(obj);
            s.epsilon = obj.epsilon;
        end
        
        function cu = upperBounds(obj,n)
            switch obj.parameters{3}
                case 1 %logsumexp approx
                    %Validate parameters
                    if obj.parameters{1} <= 0 && isinf(obj.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                        cu = [];
                    elseif obj.parameters{2} == Inf %Only min clusterdose
                        cu = Inf;
                    elseif obj.parameters{1} <= 0 %Only max clusterdose
                        cu = obj.parameters{2};
                    else %both are set sensible
                        cu = [Inf; obj.parameters{2}];
                    end
                case 2 %voxelwise
                    cu = obj.parameters{2}*ones(n,1);
                otherwise
                    error(['Min/max dose constraint evaluation method not known!']);
            end
            %cu = [Inf; obj.parameters{2}];
        end
        function cl = lowerBounds(obj,n)
            switch obj.parameters{3}
                case 1 %logsumexp approx
                    if obj.parameters{1} <= 0 && isinf(obj.parameters{2})
                        cl = [];
                    elseif obj.parameters{2} == Inf
                        cl = obj.parameters{1};
                    elseif obj.parameters{1} <= 0
                        cl = 0;
                    else
                        cl = [obj.parameters{1}; 0];
                    end
                case 2
                    cl = obj.parameters{1}*ones(n,1);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
        end
        
        function jStruct = getClusterDoseConstraintJacobianStructure(obj,n)
            switch obj.parameters{3}
                case 1 %logsumexp approx
                    %Validate parameters
                    if obj.parameters{1} <= 0 && isinf(obj.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                        jStruct = ones(n,0);
                    elseif obj.parameters{1} > 0 && isfinite(obj.parameters{2}) %both are set sensible
                        jStruct = ones(n,2);
                    else %Only min or max clusterdose
                        jStruct = ones(n,1);
                    end
                    %jStruct = ones(n,2);
                case 2
                    jStruct = speye(n);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
            
        end
        
        %% Calculates the Constraint Function value
        function cClusterDose = computeClusterDoseConstraintFunction(obj,clusterdose)
            %cDose(2) = dose_max + obj.epsilon * log( sum(exp((dose - dose_max)/obj.epsilon)));
            %cDose(1) = dose_min - obj.epsilon * log( sum(exp((dose_min - dose)/obj.epsilon)));
            switch obj.parameters{3}
                case 1 %logsumexp approx
                    cClusterDose = obj.computeClusterDoseConstraintFunctionLogSumExp(clusterdose);
                case 2
                    cClusterDose = obj.computeClusterDoseConstraintFunctionVoxelwise(clusterdose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
        end
        
        %% Calculates the Constraint jacobian
        function cClusterDoseJacob  = computeClusterDoseConstraintJacobian(obj,clusterdose)
            switch obj.parameters{3}
                case 1 %logsumexp approx
                    cClusterDoseJacob = obj.computeClusterDoseConstraintJacobianLogSumExp(clusterdose);
                case 2
                    cClusterDoseJacob = obj.computeClusterDoseConstraintJacobianVoxelwise(clusterdose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Min/max dose constraint evaluation method not known!');
            end
        end
    end
    
    methods (Access = private)
        % LogSumExp Approximation
        function cClusterDose = computeClusterDoseConstraintFunctionLogSumExp(obj,clusterdose)
            clusterdose_min = min(clusterdose);
            clusterdose_max = max(clusterdose);
            
            %Validate parameters
            if obj.parameters{1} <= 0 && isinf(obj.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                cClusterDose = [];
            elseif obj.parameters{2} == Inf %Only min dose
                cClusterDose = clusterdose_min - obj.epsilon * log( sum(exp((clusterdose_min - clusterdose)/obj.epsilon)));
            elseif obj.parameters{1} <= 0 %Only max dose
                cClusterDose = clusterdose_max + obj.epsilon * log( sum(exp((clusterdose - clusterdose_max)/obj.epsilon)));
            else %both are set sensible
                cClusterDose(2,1) = clusterdose_max + obj.epsilon * log( sum(exp((clusterdose - clusterdose_max)/obj.epsilon)));
                cClusterDose(1,1) = clusterdose_min - obj.epsilon * log( sum(exp((clusterdose_min - clusterdose)/obj.epsilon)));
            end
            
        end
        function cClusterDoseJacob  = computeClusterDoseConstraintJacobianLogSumExp(obj,clusterdose)
            %Validate parameters
            if obj.parameters{1} <= 0 && isinf(obj.parameters{2}) %Constraint doesn't make sense (min = 0 & max = Inf)
                cClusterDoseJacob = [];
            elseif obj.parameters{2} == Inf %Only min dose
                cClusterDoseJacob(:,1) = exp( (min(clusterdose)-clusterdose)/obj.epsilon );
                cClusterDoseJacob(:,1) = cClusterDoseJacob(:,1)/sum(cClusterDoseJacob(:,1));
            elseif obj.parameters{1} <= 0 %Only max dose
                cClusterDoseJacob(:,1) = exp( (clusterdose-max(clusterdose))/obj.epsilon );
                cClusterDoseJacob(:,1) = cClusterDoseJacob(:,1)/sum(cClusterDoseJacob(:,1));
            else %both are set sensible
                cClusterDoseJacob(:,1) = exp( (min(clusterdose)-clusterdose)/obj.epsilon );
                cClusterDoseJacob(:,1) = cClusterDoseJacob(:,1)/sum(cClusterDoseJacob(:,1));
                
                cClusterDoseJacob(:,2) = exp( (clusterdose-max(clusterdose))/obj.epsilon );
                cClusterDoseJacob(:,2) = cClusterDoseJacob(:,2)/sum(cClusterDoseJacob(:,2));
            end
            
            
        end
        
        %Exact voxel-wise
        function cClusterDose = computeClusterDoseConstraintFunctionVoxelwise(obj,clusterdose)
            cClusterDose = clusterdose;
        end
        function cClusterDoseJacob  = computeClusterDoseConstraintJacobianVoxelwise(obj,clusterdose)
            cClusterDoseJacob = speye(numel(clusterdose),numel(clusterdose));
        end
    end
    
end


