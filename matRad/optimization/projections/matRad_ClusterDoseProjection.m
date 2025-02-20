classdef matRad_ClusterDoseProjection < matRad_BackProjection
% matRad_DoseProjection class to compute physical dose during optimization
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    methods
        function obj = matRad_ClusterDoseProjection()
            
        end
    end
    
    methods 
        function d = computeSingleScenario(~,dij,scen,w)
            if ~isempty(dij.mClusterDose{scen})
                d = dij.mClusterDose{scen}*w;
            else
                d = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end 
        end
        
        function [dExp,dOmegaV] = computeSingleScenarioProb(~,dij,scen,w)
            if ~isempty(dij.mClusterDoseExp{scen})
                dExp = dij.mClusterDoseExp{scen}*w;
                
                for i = 1:size(dij.physicalDoseOmega,2)
                   dOmegaV{scen,i} = dij.mClusterDoseOmega{scen,i} * w;
                end 
            else
                dExp = [];
                dOmegaV = [];
            end             
        end
        
        function wGrad = projectSingleScenarioGradient(~,dij,clusterDoseGrad,scen,~)
            if ~isempty(dij.mClusterDose{scen})
                wGrad = (clusterDoseGrad{scen}' * dij.mClusterDose{scen})';
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
        
        function wGrad = projectSingleScenarioGradientProb(~,dij,dExpGrad,dOmegaVgrad,scen,~)
            if ~isempty(dij.mClusterDoseExp{scen})
                wGrad = (dExpGrad{scen}' * dij.mClusterDoseExp{scen})';
                wGrad = wGrad + 2 * dOmegaVgrad;
            else
                wGrad = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end
    end
end

