classdef (Abstract) matRad_DistributionQuantity < matRad_OptimizationQuantity

    properties
        useScenarios;
    end

    methods
        function this = matRad_DistributionQuantity(dij)
            this@matRad_OptimizationQuantity();
            if exist('dij', 'var')
                this.initializeProperties(dij);
            end

        end

        function quantityOutput = getResult(this,dij,w)
            if ~isequal(this.wCache,w)
                this.d(this.useScenarios) = arrayfun(@(scen) this.computeQuantity(dij,scen,w), this.useScenarios, 'UniformOutput',false);
                this.wCache = w;
            end
            quantityOutput = this.d;
        end

        function gradOutput = getProjectedGradient(this,dij,fGrad,w)
            if ~isequal(this.wGradCache,w)
                this.wGrad(this.useScenarios) = arrayfun(@(scen) this.projectGradient(dij,scen,fGrad,w), this.useScenarios, 'UniformOutput',false);
                this.wGradCache = w;
            end
            gradOutput = this.wGrad;
        end

        function constJacobianOutput = getProjectedJacobian(this,dij,fJacob,w)
            if ~isequal(this.wConstJacobianCache,w)
                this.wJacob{1} = this.projectConstraintJacobian(dij,fJacob,w);
                this.wConstJacobianCache = w;
            end
            constJacobianOutput = this.wJacob;
        end

        function initializeProperties(this,dij)
            
            % This is quite a mess
            if isfield(dij, 'physicalDose') && ~isempty(dij.physicalDose{1})
                distributionQuantity = 'physicalDose';
            elseif isfield(dij, 'physicalDoseExp') && ~isempty(dij.physicalDoseExp{1})
                distributionQuantity = 'physicalDoseExp';
            elseif isfield(dij, 'mAlphaDose') && ~isempty(dij.mAlphaDose{1})
                distributionQuantity = 'mAlphaDose';
            elseif isfield(dij, 'mAlphaDoseExp') && ~isempty(dij.mAlphaDoseExp{1})
                distributionQuantity = 'mAlphaDoseExp';
            end

            this.d     = cell(size(dij.(distributionQuantity)));
            this.wGrad = cell(size(dij.(distributionQuantity)));
            this.wJacob = cell(1);
        end

    end

    % methods (Abstract)
    %     function quantityOutput = computeQuantity(this, dij, scen,w)
    %         matRad_cfg = MatRad_Config.instance();
    %         matRad_cfg.dispError('The optimization quantity must implement the ''computeQuantity'' method.');
    %     end
    % 
    % 
    %     function gradientOutput = projectGradient(~,dij,scen,fGrad,w)
    %         matRad_cfg = MatRad_Config.instance();
    %         matRad_cfg.dispError('The optimization quantity must implement the ''projectGradient'' method.');
    %     end
    % 
    %     function constJacobianOutput = projectConstraintJacobian(this,dij,fJacob,w)
    %         % matRad_cfg = MatRad_Config.instance();
    %         % matRad_cfg.dispError('The optimization quantity must implement the ''projectConstraintJacobian'' method.');
    % 
    %     end
    % end
end