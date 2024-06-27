classdef matRad_cdSquaredOverdosing < ClusterDoseObjectives.matRad_ClusterDoseObjective
% matRad_SquaredOverdosing Implements a penalized squared overdosing objective
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
        name = 'Squared Overdosing';
        parameterNames = {'cd^{max}'};
        parameterTypes = {'clusterdose'};
    end
    
    properties
        parameters = {30};
        penalty = 1;
    end
    
    methods
        function obj = matRad_cdSquaredOverdosing(penalty,cdMax)
            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            %Call Superclass Constructor (for struct initialization)
            obj@ClusterDoseObjectives.matRad_ClusterDoseObjective(inputStruct);
            
            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin == 2 && isscalar(cdMax)
                    obj.parameters{1} = cdMax;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
        
        %% Calculates the Objective Function value
        function fDose = computeClusterDoseObjectiveFunction(obj,clusterDose)
            % overdose : dose minus prefered dose
            % overdose = dose - obj.parameters{1};
            overdose = (clusterDose - obj.parameters{1});
            
            % apply positive operator
            overdose(overdose<0) = 0;
            
            % claculate objective function
            fDose = obj.penalty/numel(clusterDose) * (overdose'*overdose);
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeClusterDoseObjectiveGradient(obj,clusterDose)
            % overdose : dose minus prefered dose
            % overdose = dose - obj.parameters{1};
            overdose = (clusterDose - obj.parameters{1});
            
            % apply positive operator
            overdose(overdose<0) = 0;
            
            % calculate delta
            fDoseGrad = 2 * obj.penalty/numel(clusterDose) * overdose;
        end
    end
    
end
