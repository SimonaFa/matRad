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
    end

end