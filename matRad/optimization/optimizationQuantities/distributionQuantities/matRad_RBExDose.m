classdef matRad_RBExDose < matRad_DistributionQuantity

    properties (Constant)
        quantityName = 'RBExDose';
        requiredSubquantities = {'effect'};

    end

    methods
        function this = matRad_RBExDose(dij)
            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(this, dij, scen,w)
            
            effectQt = this.getSubQuantity('effect');
            effect = effectQt.getResult(dij,w);

            quantityOutput = zeros(dij.doseGrid.numOfVoxels,1);
            [ctScen,~] = ind2sub(size(dij.physicalDose),scen);
            quantityOutput(dij.ixDose{ctScen}) = sqrt((effect{scen}(dij.ixDose{ctScen})./dij.bx{ctScen}(dij.ixDose{ctScen}))+(dij.gamma{ctScen}(dij.ixDose{ctScen}).^2)) - dij.gamma{ctScen}(dij.ixDose{ctScen});
            
        end

        function gradientOutput = projectGradient(this,dij,scen,fGrad,w)
            
            effectQt = this.getSubQuantity('effect');

            currRBEvalue = this.getResult(dij,w);

            [ctScen,~] = ind2sub(size(dij.physicalDose),scen);
 
            scaledEffect = currRBEvalue{scen} + dij.gamma{ctScen};

            fGradTemp = zeros(dij.doseGrid.numOfVoxels,1);
            fGradTemp(dij.ixDose{ctScen}) = fGrad{scen}(dij.ixDose{ctScen}) ./ (2*dij.bx{ctScen}(dij.ixDose{ctScen}).*scaledEffect(dij.ixDose{ctScen}));

            fGradRBE = cell(numel(this.useScenarios),1);
            fGradRBE{scen} = fGradTemp; 
            
            gradientOutput = effectQt.projectGradient(dij,scen,fGradRBE,w);
        end

        function constJacobianOutput = projectConstraintJacobian(this,dij,fJacob,w)
            
            effectSubQuantity = this.getSubQuantity('effect');

            currRBEvalue = this.getResult(dij,w);

            [ctScen,~] = ind2sub(size(dij.physicalDose),1);
 
            scaledEffect = currRBEvalue{1} + dij.gamma{ctScen};

            nConst = size(fJacob{1},2);
            fJacobTemp = sparse(zeros(dij.doseGrid.numOfVoxels,nConst));
            fJacobTemp(dij.ixDose{ctScen},:) = fJacob{1}(dij.ixDose{ctScen},:) ./ (2*dij.bx{ctScen}(dij.ixDose{ctScen}).*scaledEffect(dij.ixDose{ctScen}));

            constJacobianOutput = effectSubQuantity.projectConstraintJacobian(dij,{fJacobTemp},w);
        end
    end
end