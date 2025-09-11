classdef matRad_TabulatedQuantityModel < matRad_LQBasedModel

    properties
        includedFragments;
    end

    methods
        function this = matRad_TabulatedQuantityModel
            this@matRad_LQBasedModel()
        end

        function interpolatedTable = interpolateQuantityOnTables(this,referenceTable, tableToInterpolate, quantities)

            interpolatedTable = [];
            for fragmentIdx = 1:numel(tableToInterpolate)

                Z = tableToInterpolate(fragmentIdx).Z;
                A = tableToInterpolate(fragmentIdx).A;

                currReferenceTable       = this.extractFragmentsWithZA(Z, A, referenceTable);
                currTableToInterpolate   = this.extractFragmentsWithZA(Z, A, tableToInterpolate);

                for quantity = quantities
                    % For now we assume that reference table is spectra
                    % table with an "energyBin" field. This need to be made
                    % consistent later on.
                    currTable.(quantity{1}) = interp1(currTableToInterpolate.energies, currTableToInterpolate.(quantity{1}), currReferenceTable.energyBin','linear','extrap');
                end
                currTable.energies = currReferenceTable.energyBin;
                
                interpolatedTable = [interpolatedTable, currTable];
            end
        end

        function vTissueIndex = getTissueInformation(this,machine, cstDownsampled, dij,vAlphaX, ~, VdoseGrid, VdoseGridScenIx)

            matRad_cfg = MatRad_Config.instance();

            numOfCtScen = numel(vAlphaX);
            
            tmpScenVdoseGrid = cell(numOfCtScen,1);

            for s = 1:numOfCtScen            
                tmpScenVdoseGrid{s} = VdoseGrid(VdoseGridScenIx{s});
                vTissueIndex{s}     = zeros(size(tmpScenVdoseGrid{s},1),1);
            end

            if isfield(machine.data,'alphaX') && isfield(machine.data,'betaX')
                for i = 1:size(cstDownsampled,1)
                    % check if cst is compatiable
                    if ~isempty(cstDownsampled{i,5}) && isfield(cstDownsampled{i,5},'alphaX') && isfield(cstDownsampled{i,5},'betaX')
     
                        % check if base data contains alphaX and betaX
                        % (should already be checked before to instantiate the model class)
                        IdxTissue = find(ismember(machine.data(1).alphaX,cstDownsampled{i,5}.alphaX) & ...
                            ismember(machine.data(1).betaX,cstDownsampled{i,5}.betaX));
     
                        % check consitency of biological baseData and cst settings
                        if ~isempty(IdxTissue)

                            for s = 1:numOfCtScen
                                tmpScenVdoseGrid = VdoseGrid(VdoseGridScenIx{s});
                                isInVdoseGrid = ismember(tmpScenVdoseGrid,cstDownsampled{i,4}{s});
                                vTissueIndex{s}(isInVdoseGrid) = IdxTissue;
                            end
                        else
                            matRad_cfg.dispError('Biological base data and cst are inconsistent!');
                        end
     
                    else
                        for s = 1:numOfCtScen
                            vTissueIndex{s}(:) = 1;
                        end
                        matRad_cfg.dispWarning('\tTissue type of %s was set to 1\n',cstDownsampled{i,2});
                    end
                end
                
                matRad_cfg.dispInfo('done.\n');
             else
                matRad_cfg.dispError('Base data is missing alphaX and/or betaX!');
            end

        end

    end


    methods (Static)

        % function filteredTable = extractFragmentsFromTable(table, fragments)
        % 
        %     matRad_cfg = MatRad_Config.instance();
        % 
        %     % Check if all fragments were selected from the table
        %     if ~all(ismember([[fragments.Z]' [fragments.A]'], [[filteredTable.Z]' [filteredTable.A]'], 'rows'))
        %         matRad_cfg.dispError(sprintf('At least one of the required fragments from table:%s is not available', table.meta.name));
        %     end
        % end


        function outData = extractFragmentsWithZA(Z,A, structArray)

            availableZ = [structArray.Z]';
            availableA = [structArray.A]';

            availableFragments = [availableZ, availableA];
            
            indexes = ismember(availableFragments,[Z',A'], 'rows');

            outData = structArray(indexes);
            
        end
    end
end