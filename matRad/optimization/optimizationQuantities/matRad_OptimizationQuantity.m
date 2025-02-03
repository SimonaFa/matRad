classdef (Abstract) matRad_OptimizationQuantity < handle

    properties (Abstract,Constant)
        quantityName;
        requiredSubquantities;
    end

    properties
        subQuantities;
        d;
        wGrad;
        wCache;
        wGradCache;
        %c;
        wJacob;
        %wConstraintCache;
        wConstJacobianCache;
    end

    methods
        function this = matRad_OptimizationQuantity()

        end

        % Implemented by distibution/scalar quantities
        function output = getResult(~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Function needs to be implemented by subclass');
            output = [];

        end

        function output = getProjectedGradient(~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Function needs to be implemented by subclass');
            output = [];
        end

        % function output = getConstraintResult(~)
        %     matRad_cfg = MatRad_Config.instance();
        %     matRad_cfg.dispError('Function needs to be implemented by subclass');
        %     output = [];
        % end

        function output = getProjectedJacobian(~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Function needs to be implemented by subclass');
            output = [];
        end

        % Implemented by specific subclass
        function output = computeQuantity(~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Function needs to be implemented by subclass');
            output = [];
        end

         function output = projectGradient(~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Function needs to be implemented by subclass');
            output = [];
         end

         % function output = computeConstraint(~)
         %    matRad_cfg = MatRad_Config.instance();
         %    matRad_cfg.dispError('Function needs to be implemented by subclass');
         %    output = [];
         % end

         function output = projectConstraintJacobian(~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Function needs to be implemented by subclass');
            output = [];
         end



         function subQuantityInstance = getSubQuantity(this, name)
            subQuantityInstance = this.subQuantities{cellfun(@(x) strcmp(x.quantityName, name), this.subQuantities)};
         end
    end

    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc,alphaX,betaX)

        end
    end

end