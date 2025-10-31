classdef matRad_Fluence < matRad_DijDistributionQuantity

    properties (Constant)
        quantityName = 'Fluence';
        requiredSubquantities = {};
        
        dijField = {'mFluence'};
    end

    methods

        function this = matRad_Fluence(dij)

            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DijDistributionQuantity(supArg{:});
        end
    end

end