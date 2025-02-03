classdef matRad_Effect < matRad_DistributionQuantity

    properties (Constant)
        quantityName = 'effect';
        requiredSubquantities = {'AlphaDose', 'SqrtBetaDose'};
    end

    properties
        
    end

    methods
        function this = matRad_Effect(dij)
            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});

        end

        function quantityOutput = computeQuantity(this, dij, scen,w)

            alphaDoseQt = this.getSubQuantity('AlphaDose');
            alphaDose = alphaDoseQt.getResult(dij,w);

            betaDoseQt = this.getSubQuantity('SqrtBetaDose');
            betaDose   = betaDoseQt.getResult(dij,w);

            quantityOutput = alphaDose{scen} + betaDose{scen}.^2;
        end

        function gradientProjectionOutput = projectGradient(this,dij,scen,fGrad,w)
            
            % These functions call directly the subquantity calculation,
            % the w cache is not checked by the subquantity and the result
            % is not stored there, it is only stored by this quantity if
            % needed.
            alphaDoseQt = this.getSubQuantity('AlphaDose');
            betaDoseQt = this.getSubQuantity('SqrtBetaDose');

            alphaGrad = alphaDoseQt.projectGradient(dij,scen,fGrad,w);

            % Get result is instead passed through the w cache check of the
            % subfunction, it is a quantity, so depends only on the current
            % status of the weigts. It does not make sense to have
            % different quantities being called with different weights in
            % the same optimization step, so the heavy calculation here
            % should only be executed on the first call to this function
            % within the optimmization step
            betaDose  = betaDoseQt.getResult(dij,w);
            
            fBetaGrad = cell(numel(this.useScenarios),1);
            fBetaGrad{scen} = 2 * (fGrad{scen} .* betaDose{scen});
            
            betaGrad = betaDoseQt.projectGradient(dij,scen,fBetaGrad,w);

            vBias = alphaGrad;
            mPsi  = betaGrad;

            gradientProjectionOutput = (vBias + mPsi);

        end

        function constJacobianOutput = projectConstraintJacobian(this,dij,fJacob,w)
            
            alphaSubQuantity    = this.getSubQuantity('AlphaDose');
            sqrtBetaSubQuantity = this.getSubQuantity('SqrtBetaDose');
            
            alphaJacob    = alphaSubQuantity.projectConstraintJacobian(dij,fJacob);
            %sqrtBetaJacob = sqrtBetaSubQuantity.projectConstraintJacobian(dij,fJacob);
            sqrtBetaDose  = sqrtBetaSubQuantity.getResult(dij,w);

            fJacobBeta = {2 * sqrtBetaDose{1} .* fJacob{1}};
            betaJacob = sqrtBetaSubQuantity.projectConstraintJacobian(dij,fJacobBeta);
            
            constJacobianOutput = alphaJacob + betaJacob;
        end
    
    end

    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc,alphaX,betaX)
            doses = optiFunc.getDoseParameters();
            effect = alphaX*doses + betaX*doses.^2;
            optiFunc = optiFunc.setDoseParameters(effect);
        end
    end

end