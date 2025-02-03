classdef matRad_LETd < matRad_DijDistributionQuantity

    properties (Constant)
        quantityName = 'LETd';
        requiredSubquantities = {};
        
        dijField = {'mLETDose'};
    end

    methods

        function this = matRad_LETd(dij)

            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DijDistributionQuantity(supArg{:});
        end
    end

end