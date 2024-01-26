classdef matRad_cdSquaredDeviation < ClusterDoseObjectives.matRad_ClusterDoseObjective
% matRad_SquaredDeviation Implements a penalized least squares objective
%   See matRad_DoseObjective for interface description
%
% References 
%     -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
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
        name = 'Squared Deviation';
        parameterNames = {'cd^{ref}'};
        parameterTypes = {'clusterdose'};
    end
    
    properties
        parameters = {10};
        penalty = 1;
    end
    
    methods
        function obj = matRad_cdSquaredDeviation(penalty,dRef)
            
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
                if nargin == 2 && isscalar(dRef)
                    obj.parameters{1} = dRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
   
        %% Calculates the Objective Function value
        function fClusterDose = computeClusterDoseObjectiveFunction(obj,clusterDose)
            % deviation : dose minus prefered dose
            %deviation = dose - obj.parameters{1};
            deviation = (clusterDose - obj.parameters{1});
            
            % claculate objective function
            fClusterDose = obj.penalty/numel(clusterDose) * (deviation'*deviation);
        end
        
        %% Calculates the Objective Function gradient
        function fClusterDoseGrad   = computeClusterDoseObjectiveGradient(obj,clusterDose)
            % deviation : Dose minus prefered dose
            % deviation = dose - obj.parameters{1};
            deviation = (clusterDose - obj.parameters{1});
            
            % calculate delta
            fClusterDoseGrad = 2 * obj.penalty/numel(clusterDose) * deviation;
        end
    end
    
end
