classdef (Abstract) matRad_DijDistributionQuantity < matRad_DistributionQuantity

    properties (Abstract, Constant)
        dijField;
    end

    methods
        function this = matRad_DijDistributionQuantity(dij)
            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DistributionQuantity(supArg{:});

            if nargin>0 && ~isfield(dij.(this.dijField{1}))

                matrad_cfg = MatRad_Config().instance();
                matrad_cfg.dispError(['Provided dij does not contain required field name: ', dij.dijField{1}]);
            end
        
        end

        function quantityOutput = computeQuantity(this, dij, scen,w)
        % This function implements the product dij * w for a generic
        % field of the dij matrix to be instantiated by the specific
        % subclass

            % Check that the field is present and not empty for the
            % required scenario
            if ~isempty(dij.(this.dijField{1}){scen})
        
                % Compute the product
                quantityOutput = dij.(this.dijField{1}){scen}*w;
        
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
                quantityOutput = [];
            end
        end

        function gradientOutput = projectGradient(this,dij,scen,fGrad,w)
        % This function implements the product of the objective function 
        % gradient and the quantity specific gradient for a generic
        % field of the dij matrix to be instantiated by the specific
        % subclass
    
            if ~isempty(dij.(this.dijField{1}){scen})
                gradientOutput = (fGrad{scen}' * dij.(this.dijField{1}){scen})';
            else
                gradientOutput = [];
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Empty scenario in optimization detected! This should not happen...\n');
            end
        end

        function constJacobianOutput = projectConstraintJacobian(this,dij,fJacob,~)
        % This function implements the product of the constraint function 
        % jacobian and the quantity specific gradient for a generic
        % field of the dij matrix to be instantiated by the specific
        % subclass
            constJacobianOutput = fJacob{1}' * dij.(this.dijField{1}){1};
        end


    end

end