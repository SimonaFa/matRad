classdef matRad_MeanClusterDose < DoseObjectives.matRad_ClusterDoseObjective
% matRad_MeanDose Implements a penalized MeanDose objective
%   See matRad_DoseObjective for interface description
%
% References
%   -
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
        name = 'Mean Cluster Dose';
        parameterNames = {'cd^{ref}','f_{diff}'}; %When optimizing to a reference, one might consider using a quadratic relationship with a non-linear optimizer
        parameterTypes = {'clusterdose',{'Linear','Quadratic'}};
    end
    
    properties
        parameters = {0,1};        
        penalty = 1;
    end
    
    methods 
        function obj = matRad_MeanClusterDose(penalty,cdMeanRef,fDiff)
           
            % if we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@DoseObjectives.matRad_ClusterDoseObjective(inputStruct);
            
            if ~initFromStruct
                if nargin < 3 || ~ischar(fDiff)
                    fDiff = 'Linear';
                end
                
                fDiffIx = find(strcmp(fDiff,obj.parameterTypes{2}));
                
                if isempty(fDiffIx) || numel(fDiffIx) > 1
                    fDiffIx = 1;
                    matRad_cfg = MatRad_Config.instance();                    
                    matRad_cfg.dispWarning('Mean dose difference function can only be %s! Using %s difference.', strjoin(obj.parameterTypes{2},' or '), obj.parameterTypes{2}{fDiffIx});
                end
                
                obj.parameters{2} = fDiffIx;


                if nargin >= 2 && isscalar(cdMeanRef)
                    obj.parameters{1} = cdMeanRef;
                end

                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end

            %% Downwards compatibility / set default values
            %TODO: maybe move into set method for parameters
            if numel(obj.parameters) < 1
                obj.parameters{1} = 0;
            end

            if numel(obj.parameters) < 2
                obj.parameters{2} = 1;
            end
            
        end       
        
        %% Calculates the Objective Function value
        function fDose = computeClusterDoseObjectiveFunction(obj,dose)
            switch obj.parameters{2}
                case 1
                    fDose = obj.penalty * obj.objectiveLinearDiff(dose);
                case 2
                    fDose = obj.penalty * obj.objectiveQuadraticDiff(dose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean Dose Objective!',obj.parameterNames{2});  
            end
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeClusterDoseObjectiveGradient(obj,clusterDose)
            switch obj.parameters{2}
                case 1
                    fDoseGrad = obj.penalty * obj.gradientLinearDiff(clusterDose);
                case 2
                    fDoseGrad = obj.penalty * obj.gradientQuadraticDiff(clusterDose);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid setting for %s in Mean Dose Objective!',obj.parameterNames{2});  
            end
        end
    end

    methods (Access = protected)
        function fDose = objectiveQuadraticDiff(obj,clusterDose)
            fDose = (mean(clusterDose(:)) - obj.parameters{1})^2;
        end

        function fDoseGrad = gradientQuadraticDiff(obj,clusterDose)
            fDoseGrad = 2*(mean(clusterDose(:))-obj.parameters{1}) * ones(size(clusterDose(:)))/numel(clusterDose);
        end

        function fDose = objectiveLinearDiff(obj,clusterDose)
            fDose = abs(mean(clusterDose(:)) - obj.parameters{1});
        end

        function fDoseGrad = gradientLinearDiff(obj,clusterDose)
            fDoseGrad = (1/numel(clusterDose))*sign(clusterDose(:)-obj.parameters{1});
        end
    end
    
end

