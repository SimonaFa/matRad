classdef matRad_BED < matRad_DistributionQuantity

    properties (Constant)
        quantityName = 'BED';
        requiredSubquantities = {'effect'};
    end

    methods
        function this = matRad_BED(dij)
            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(this, dij, scen,w)

            [ctScen,~] = ind2sub(size(dij.physicalDose),scen);

            effectSubquantity = this.getSubQuantity('effect');
            effect = effectSubquantity.getResult(dij,w);
            
            quantityOutput = zeros(dij.doseGrid.numOfVoxels,1);
            quantityOutput(dij.ixDose{ctScen}) = effect{scen}(dij.ixDose{ctScen})./ dij.ax{ctScen}(dij.ixDose{ctScen});
 
        end

        function gradientProjectionOutput = projectGradient(this,dij,scen,fGrad,w)
            
            [ctScen,~] = ind2sub(size(dij.physicalDose),scen);

            fGradtmp{scen} = zeros(size(fGrad{scen}));
            fGradtmp{scen}(dij.ixDose{ctScen}) = fGrad{scen}(dij.ixDose{ctScen})./dij.ax{scen}(dij.ixDose{ctScen});
            
            effectSubquantity = this.getSubQuantity('effect');

            gradientProjectionOutput = effectSubquantity.projectGradient(dij,scen,fGradtmp,w);
        end
    end

    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc,alphaX,betaX)
            doses = optiFunc.getDoseParameters();    
            BED = doses*(1 + doses/(alphaX/betaX));        
            optiFunc = optiFunc.setDoseParameters(BED);
        end
    end
end