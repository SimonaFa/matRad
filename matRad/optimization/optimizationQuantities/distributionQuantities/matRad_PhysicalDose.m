classdef matRad_PhysicalDose < matRad_DijDistributionQuantity

    properties (Constant)
        quantityName = 'physicalDose';
        requiredSubquantities = {};
        
        dijField = {'physicalDose'};

        doseFields = {'physicalDose','doseToWater'};
        doseQuantities = {'','_std','_batchStd'}; 
    end

    methods

        function this = matRad_PhysicalDose(dij)

            if nargin>0
                supArg = {dij};
            else
                supArg = {};
            end
            
            this@matRad_DijDistributionQuantity(supArg{:});
        end

        function resultGUI = matRad_calcQuantityCubes(~, dij, resultGUI, scenNum, beamInfo)
            % compute physical dose for all beams individually and together
            for j = 1:length(matRad_PhysicalDose.doseFields)
                for k = 1:length(matRad_PhysicalDose.doseQuantities)
                    % Check if combination is a field in dij, otherwise skip
                    if isfield(dij,[matRad_PhysicalDose.doseFields{j} matRad_PhysicalDose.doseQuantities{k}])
                        % Handle standard deviation fields and add quadratically
                        if ~isempty(strfind(lower(matRad_PhysicalDose.doseQuantities{1}),'std'))
                            for i = 1:length(beamInfo)
                                resultGUI.([matRad_PhysicalDose.doseFields{j}, matRad_PhysicalDose.doseQuantities{k}, beamInfo(i).suffix]) = sqrt(reshape(full(dij.([matRad_PhysicalDose.doseFields{j} matRad_PhysicalDose.doseQuantities{k}]){scenNum}.^2 * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions));
                                resultGUI.([matRad_PhysicalDose.doseFields{j}, matRad_PhysicalDose.doseQuantities{k}, beamInfo(i).suffix])(isnan(resultGUI.([matRad_PhysicalDose.doseFields{j}, matRad_PhysicalDose.doseQuantities{k}, beamInfo(i).suffix]))) = 0;
                            end
                            % Handle normal fields as usual
                        else
                            for i = 1:length(beamInfo)
                                resultGUI.([matRad_PhysicalDose.doseFields{j}, matRad_PhysicalDose.doseQuantities{k}, beamInfo(i).suffix]) = reshape(full(dij.([matRad_PhysicalDose.doseFields{j} matRad_PhysicalDose.doseQuantities{k}]){scenNum} * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions);
                            end
                        end
                    end
                end
            end
        end

    end

end