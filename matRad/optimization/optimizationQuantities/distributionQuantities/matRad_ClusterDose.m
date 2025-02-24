classdef matRad_ClusterDose < matRad_DijDistributionQuantity

    properties (Constant)
        quantityName = 'clusterDose';
        requiredSubquantities = {};
        
        dijField = {'mClusterDose'};
    end

    methods

        function this = matRad_ClusterDose(dij)

            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DijDistributionQuantity(supArg{:});
        end
%{

        function resultGUI = matRad_calcQuantityCubes(~, dij, resultGUI, scenNum, beamInfo)
            if isfield(dij,'mClusterDose')
                for i = 1:length(beamInfo)
                    resultGUI.(['clusterDose', beamInfo(i).suffix])     = reshape(full(dij.mClusterDose{scenNum} * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions);
                end
            end
        end
%}
    end

end