classdef matRad_SqrtBetaDose < matRad_DijDistributionQuantity

    properties (Constant)
        quantityName = 'SqrtBetaDose';
        requiredSubquantities = {};

        dijField = {'mSqrtBetaDose'};
    end

    properties
        parameter;
    end
    methods
        function this = matRad_SqrtBetaDose(dij)
            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DijDistributionQuantity(supArg{:});
        end
    end
end