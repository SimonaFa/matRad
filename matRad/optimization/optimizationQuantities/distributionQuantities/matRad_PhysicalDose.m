classdef matRad_PhysicalDose < matRad_DijDistributionQuantity

    properties (Constant)
        quantityName = 'physicalDose';
        requiredSubquantities = {};
        
        dijField = {'physicalDose'};
    end

    methods

        function this = matRad_PhysicalDose(dij)

            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DijDistributionQuantity(supArg{:});
        end
    end

end