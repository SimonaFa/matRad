classdef (Abstract) matRad_TabulatedDoseAveragedKernelModel < matRad_TabulatedQuantityModel

    properties (Abstract, Constant)
        quantitiesToAverage;
    end
    properties
        ionsName = {'H1', 'He', 'Li', 'Be', 'B', 'C'};
        defaultZ    = [1, 2, 3, 4, 5, 6];
        defaultA   = [1, 4, 7, 9, 11, 12];
    end

    properties
        stoppingPowerTable;
        stoppingPowerTableName;
        quantityTable;
        quantityTableName;
    end

    methods
        function this = matRad_TabulatedDoseAveragedKernelModel
            this@matRad_TabulatedQuantityModel()
        end

        function machine = computeKernels(this, machine)
            % This is the main function.

            % Load the stopping power table for the selected fragments
            spTable = this.getSpTable(this.includedFragments);

            % Check table consistency. TODO: move this form here
            for eIdx=1:numel(machine.data)
                for fragmentIdx=1:numel(machine.data(eIdx).Fluence.spectra)
                    if isnan(machine.data(eIdx).Fluence.spectra(fragmentIdx).A)
    
                        matRad_cfg = MatRad_Config.instance();
    
                        [~, tmpIdx] = intersect([this.includedFragments.Z], machine.data(eIdx).Fluence.spectra(fragmentIdx).Z);
    
                        if ~isempty(tmpIdx)
                            machine.data(eIdx).Fluence.spectra(fragmentIdx).A = this.includedFragments(tmpIdx).A;
                            if eIdx == 1
                                matRad_cfg.dispWarning(sprintf('NaN value for A detected for ion with Z:%d, assuming default value of A:%d', machine.data(eIdx).Fluence.spectra(fragmentIdx).Z, machine.data(eIdx).Fluence.spectra(fragmentIdx).A));
                            end
                        end
    
                    end
                end
            end



            for eIdx=1:numel(machine.data)

                currMachineData = machine.data(eIdx);

                spectra = this.getSpectraFromMachine(currMachineData, this.includedFragments, 'Fluence');

                % For now energies in the base data are in total energy,
                % move it to energy per nucleon
                % for i=1:numel(spectra)
                %     spectra(i).energyBin = spectra(i).energyBin./spectra(i).A;
                % end
                
                currentSpTable   = this.interpolateQuantityOnTables(spectra, spTable, {'dEdx'});  % Syntax: interpolateQuantityOnTables(referenceTable, tableToInterpolate, quantityToInterpolate)
                quantityToWeight = this.interpolateQuantityOnSpectra(spectra);

                kernels = this.computeDoseAveragedKernels(quantityToWeight, spectra, currentSpTable);

                for kernelName=fieldnames(kernels)'
                    machine.data(eIdx).(kernelName{1}) = kernels.(kernelName{1});
                end


            end

        end
        
        function doseAveragedKernels = computeDoseAveragedKernels(this,quantity, spectra, sp)

            
            % Get normalization for each fragment
            denominator = arrayfun(@(fragmentSpectrum, fragmentSP) fragmentSpectrum.fluenceSpectrum' * fragmentSP.dEdx, spectra, sp, 'UniformOutput',false);
            
            % Sum over the fragments
            denominator = sum([denominator{:}],2);

            % Get the numerator for every fragment and quantity
            for quantityName=this.quantitiesToAverage
                currentQuantity = quantityName{1};
                numerator.(currentQuantity) = arrayfun(@(fragmentSpectrum, fragmentSP, fragQuantity) fragmentSpectrum.fluenceSpectrum' * (fragmentSP.dEdx .* fragQuantity.(currentQuantity)), spectra, sp, quantity, 'UniformOutput',false);

                % Sum over the fragments
                numerator.(currentQuantity) = sum(cat(3,numerator.(currentQuantity){:}),3);

                % Get the kernels
                doseAveragedKernels.(currentQuantity) = numerator.(currentQuantity)./denominator;
            end


        end

         function outQuantity = interpolateQuantityOnSpectra(this, spectra)

            outQuantity = [];
            fTable = this.getTable(this.includedFragments);
            outQuantity = this.interpolateQuantityOnTables(spectra, fTable, this.quantitiesInTable);

         end

         function table = getTable(this, fragments)

            
             if isempty(this.quantityTable)
                this.quantityTable = this.loadTable(this.quantityTableName);
             end

             if ~isfield(this.quantityTable.data, 'Z') || ~isfield(this.quantityTable.data, 'A') && isfield(this.quantityTable.data, 'includedIons')
                included = this.quantityTable.data(1).includedIons;
                tmpZ = arrayfun(@(ion) this.defaultZ(strcmp(ion{1}, this.ionsName)), included);
                tmpA = arrayfun(@(ion) this.defaultA(strcmp(ion{1}, this.ionsName)), included);
                [this.quantityTable.data.Z] = deal(tmpZ);
                [this.quantityTable.data.A] = deal(tmpA);
             end

             tmpZ = [fragments.Z];
             tmpA = [fragments.A];

             if any([fragments.Z] == 1 & [fragments.A] ~= 1)
                 [tmpA([fragments.Z] == 1)] = deal(1);
             end

             table = this.extractFragmentsWithZA([tmpZ], [tmpA], this.quantityTable.data);

             if sum(([table.Z]==1))>1
                 [table.A] = deal(fragments.A);
             end

        end

        function spTable = getSpTable(this, fragments)

            if isempty(this.stoppingPowerTable)
                this.stoppingPowerTable = this.loadStoppingPowerTable(this.stoppingPowerTableName);
            end

            tmpZ = [fragments.Z];
            tmpA = [fragments.A];

            if any([fragments.Z] == 1 & [fragments.A] ~= 1)
                [tmpA([fragments.Z] == 1)] = deal(1);
            end

            spTable = this.extractFragmentsWithZA(tmpZ, tmpA, this.stoppingPowerTable.data);

            if sum(([spTable.Z]==1))>1
               [spTable.A] = deal(fragments.A);
            end

            % TODO: Check that all fragments got out
            if numel(spTable) ~= numel(fragments)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('At least one requested fragment is not available in the stopping power table');
            end
        end

        function spectra = getSpectraFromMachine(this,machineData, fragments, spectraType)
            spectra = this.extractFragmentsWithZA([fragments.Z], [fragments.A],machineData.(spectraType).spectra);
       end


    end

    methods (Static)

        function stoppingPowerTable = loadStoppingPowerTable(fileName)

            matRad_cfg = MatRad_Config.instance();

            if isempty(fileName)
                matRad_cfg.dispError('How am I supposed to compute the stopping power if you don''t give me a stopping power table!?\n Please set the stoppingPowerTableName.');
            end

            searchPath = {fullfile(matRad_cfg.matRadSrcRoot,'bioModels','StoppingPowerTables'),...    % default matrad folder
                fullfile(matRad_cfg.primaryUserFolder, 'StoppingPowerTables')};             % user defined RBE table

 
            try
                load(fullfile(searchPath{1}, [fileName, '.mat']), 'SPtable');
            catch
                try
                    load(fullfile(searchPath{2}, [fileName, '.mat']), 'SPtable');

                catch
                    matRad_cfg.dispError('Cannot find Stopping Power Table: %s', fileName);
                end
            end

            stoppingPowerTable = SPtable;
        end

    end
end