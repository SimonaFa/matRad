classdef matRad_IpSquaredxFluence < matRad_DijDistributionQuantity

    properties (Constant)
        quantityName = 'SquaredIp';
        requiredSubquantities = {'Fluence'};
        
        dijField = {'mSquaredIp'};
    end

    methods

        function this = matRad_IpSquaredxFluence(dij)

            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DijDistributionQuantity(supArg{:});
        end

    end

end