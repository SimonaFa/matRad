classdef matRad_ConstantRBExDose < matRad_DijDistributionQuantity

    properties (Constant)
        quantityName = 'constantRBExDose';
        requiredSubquantities = {};

        dijField = {'physicalDose'};
    end

    methods

        function this = matRad_ConstantRBExDose(dij)

            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DijDistributionQuantity(supArg{:});
        end

        function quantityOutput = computeQuantity(this, dij, scen,w)

            quantityOutput = computeQuantity@matRad_DijDistributionQuantity(this,dij,scen,w);
            quantityOutput = quantityOutput * dij.RBE;

        end

        function gradientOutput = projectGradient(this,dij,scen,fGrad,w)

            gradientOutput = projectGradient@matRad_DijDistributionQuantity(this,dij,scen,fGrad,w);
            gradientOutput = gradientOutput * dij.RBE;

        end

        function constJacobianOutput = projectConstraintJacobian(this,dij,fJacob,~)
            constJacobianOutput = projectConstraintJacobian@matRad_DijDistributionQuantity(this,dij,fJacob);
            constJacobianOutput = constJacobianOutput * dij.RBE;

        end

    end

end