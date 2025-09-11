classdef matRad_AlphaDose < matRad_DijDistributionQuantity

    properties (Constant)
        quantityName = 'AlphaDose';
        requiredSubquantities = {};

        dijField = {'mAlphaDose'};
    end

    methods
        function this = matRad_AlphaDose(dij)
            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DijDistributionQuantity(supArg{:});
        end
    end
end