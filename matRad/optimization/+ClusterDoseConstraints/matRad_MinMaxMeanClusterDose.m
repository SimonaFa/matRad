classdef matRad_MinMaxMeanClusterDose < ClusterDoseConstraints.matRad_ClusterDoseConstraint
    % matRad_MinMaxMeanDose Implements a MinMaxMeanDose constraint
    %   See matRad_DoseConstraint for interface description
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
        name = 'mean clusterdose constraint';
        parameterNames = {'\mu_cd^{min}', '\mu_cd^{max}'};
        %parameterIsDose = logical([1 1]);
        parameterTypes = {'clusterdose','clusterdose'};
    end
    
    properties
        parameters = {0,30};
    end
    
    methods
        function obj = matRad_MinMaxMeanClusterDose(minMeanCDose,maxMeanCDose)
            
            %If we have a struct in first argument
            if nargin == 1 && isstruct(minMeanCDose)
                initFromStruct = true;
                inputStruct = minMeanCDose;
            else
                initFromStruct = false;
                inputStruct = [];
            end
            
            % Call Superclass Constructor (for struct initialization)
            obj@ClusterDoseConstraints.matRad_ClusterDoseConstraint(inputStruct);
            
            % now handle initialization from other parameters
            if ~initFromStruct
                
                if nargin == 2 && isscalar(maxMeanCDose)
                    obj.parameters{2} = maxMeanCDose;
                end
                
                if nargin >= 1 && isscalar(minMeanCDose)
                    obj.parameters{1} = minMeanCDose;
                end
                
            end
        end
        
        function cu = upperBounds(obj,n)
            cu = obj.parameters{2};
        end
        function cl = lowerBounds(obj,n)
            cl = obj.parameters{1};
        end
        
        %Overloads the struct function to add constraint specific
        %parameters
        function s = struct(obj)
            s = struct@ClusterDoseConstraints.matRad_ClusterDoseConstraint(obj);
            %Nothing to do here...
        end
        
        
        %% Calculates the Constraint Function value
        function cClusterDose = computeDoseConstraintFunction(obj,clusterdose)
            cClusterDose = mean(clusterdose);
        end
        
        %% Calculates the Constraint jacobian
        function cClusterDoseJacob  = computeClusterDoseConstraintJacobian(obj,clusterdose)
            cClusterDoseJacob = ones(numel(clusterdose),1)./numel(clusterdose);
        end
    end
    
end


