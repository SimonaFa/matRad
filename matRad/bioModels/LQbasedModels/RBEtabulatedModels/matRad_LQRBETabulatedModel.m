classdef (Abstract) matRad_LQRBETabulatedModel < matRad_LQBasedModel
% This is an Abstract class implementig a tabulated RBE model.
% The model can handle multiple tissue alphaX/betaX ratio specified by the 
% cst structure, as long as a compatible RBEtable is provided.
%
% Properties of the model that can be defined by the user in pln.propDoseCalc.bioProperties:
%   RBEtable:       name of the specific table to be included
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    properties
        RBEtableName;
        includedFragments;          % Fragments to include from the RBE table. Can be "all", or a struct with filed "Z" or a struct withj fields "Z" and "A"
    end

    properties (SetAccess = protected, GetAccess = public)
        RBEtable;
        fragmentIndexesInBaseData;
        fragmentIndexesInTable;
 
        tissueAlphaX;               % array containing the alphaX values given in cst{i,5}. Dimension is (1,max(TissueClasses))
        tissueBetaX;                % array containing the betaX  values given in cst{i,5}. Dimension is (1,max(TissueClasses))
        availableAlphaInTable;
        availableBetaInTable;
        availableFragmentsInTable;
        defaultRBETable;
        defaultGenericFragments;
        defaultFragmentZ;
    end

    methods
        function this = matRad_LQRBETabulatedModel()
            %matRad_cfg = MatRad_Config.instance();
            
            this@matRad_LQBasedModel();
            
            this.assignDefaultProperties();
            
            % This just for testing
            %this.defaultRBETable = 'RBEtable_rapidLEM_Russo2011_longErange_LEMI30';
        
        end

        function [alphaE,betaE] = interpolateRBETableForBixel(this,interpEnergies, fragment, tissueClass)
            % This function interpolates the correct table for the input
            % fragment and tissue class
            
            fragmentRBEtable = this.selectDataTableForFragment(fragment,tissueClass);
            
            alphaE = matRad_interp1(fragmentRBEtable.energies, fragmentRBEtable.alpha, interpEnergies{1,fragment});
            betaE  = matRad_interp1(fragmentRBEtable.energies, fragmentRBEtable.beta,  interpEnergies{1,fragment});
           
            

            if min(fragmentRBEtable.energies(:)) > min(interpEnergies{1,fragment}(:))
                %matRad_cfg = MatRad_Config.instance();
                % Switch off for now
                %matRad_cfg.dispWarning('The energy is out of the RBEtable range. The min energy will be set to the min energy from the RBEtable.')
                alphaE(interpEnergies{1,fragment} < min(fragmentRBEtable.energies)) = fragmentRBEtable.alpha(1,1);
                betaE(interpEnergies{1,fragment} < min(fragmentRBEtable.energies))  = fragmentRBEtable.beta(1,1);
            end
            
            if max(fragmentRBEtable.energies(:)) < max(interpEnergies{1,fragment}(:))
                %matRad_cfg = MatRad_Config.instance();
                % Switch off for now
                %matRad_cfg.dispWarning('The energy is out of the RBEtable range. The max energy will be set to the max energy from the RBEtable.')
                alphaE(interpEnergies{1,fragment} > max(fragmentRBEtable.energies)) = fragmentRBEtable.alpha(end,1);
                betaE(interpEnergies{1,fragment} > max(fragmentRBEtable.energies))  = fragmentRBEtable.beta(end,1);
            end


        end

        function vTissueIndex = getTissueInformation(this,~, cst, dij,vAlphaX, vBetaX, VdoseGrid, VdoseGridScenIx)

            % The tissue index output corresponds to an array of
            % size(nVoxels,1) containing indexing to the data in the RBE
            % to identify the alpha/beta ratio
            matRad_cfg = MatRad_Config.instance();

            numOfCtScen = numel(vAlphaX);

            cstDownsampled = matRad_setOverlapPriorities(cst);
            
            tmpScenVdoseGrid = cell(numOfCtScen,1);

            for s = 1:numOfCtScen            
                tmpScenVdoseGrid{s} = VdoseGrid(VdoseGridScenIx{s});
                vTissueIndex{s}     = zeros(size(tmpScenVdoseGrid{s},1),1);
            end

            for s = 1:numOfCtScen
                % Get all unique alphaX/betaX values
                allAlphaBetaRatios = unique([vAlphaX{s}, vBetaX{s}], 'rows');
                for i=1:size(allAlphaBetaRatios,1)

                    % Check if this alpha beta ratio is available in the
                    % RBEtable
                    tableDataIndex = find(ismember([this.availableAlphaInTable, this.availableBetaInTable], allAlphaBetaRatios(i,:), 'rows'));

                    if isempty(tableDataIndex)
                        matRad_cfg.dispError('One or more structures in the cst have an Alpha/Beta ratio of: %1.2f/%1.2f but no data for this ratio was found in the RBE table.',allAlphaBetaRatios(i,1), allAlphaBetaRatios(i,2));
                    else
                        % Select only the voxels with the current alpha/beta
                        % ration and assign increasing tissue index
                        voxelsWithThisAlphaBetaRatio = ismember([vAlphaX{s}, vBetaX{s}], allAlphaBetaRatios(i,:), 'rows');
                        vTissueIndex{s}(voxelsWithThisAlphaBetaRatio) = tableDataIndex;
                    end
                end
            end

            % For now only considering the alpha/beta ratios of the first
            % scenarios. If different ct scenarios have different
            % alpha/beta ratios need to modify this
            allAlphaBetaRatios = unique([vAlphaX{1}, vBetaX{1}], 'rows');

            this.tissueAlphaX = allAlphaBetaRatios(:,1);
            this.tissueBetaX  = allAlphaBetaRatios(:,2);

            matRad_cfg.dispInfo('done.\n');
            
        end

        function fragmentRBEtable = selectDataTableForFragment(this,fragment,tissueClass)
            
            % This function selects the specific fragment and tissue class
            % entry in the RBE table.
            % The RBE table should be a struct containing the following
            % fields:
            %   .meta   struct with meta information
            %   .data   struct array with one entry for each tissue class
            % In turn, the data struct should contain the subfields:
            %   includedIons    ions included for the specific table (i.e 'H', 'He', 'C', ...);                                      
            %   energies        array containing the energies corresponding
            %                   to the specified alpha/beta tables
            %   alpha           matrix (#energies,#fragments) specifieng the alpha
            %                   value for the included fragments and
            %                   energies
            %   beta            matrix (#energies,#fragments) specifieng the alpha
            %                   value for the included fragments and
            %                   energies

            matRad_cfg = MatRad_Config.instance();

            fragmentRBEtable.alpha    = this.RBEtable.data(tissueClass).alpha(:,fragment);
            fragmentRBEtable.beta     = this.RBEtable.data(tissueClass).beta(:,fragment);
            fragmentRBEtable.energies = this.RBEtable.data(tissueClass).energies.*this.RBEtable.data(tissueClass).includedIons(fragment).A; % RBEtable should contain energy per nucleon
            
         
        end

        function updateFragmentIdxInTable(this,includedFragments)

            if isempty(this.RBEtable)
                return;
            end

            matRad_cfg = MatRad_Config.instance();

            if ischar(includedFragments)
                if strcmp(includedFragments, 'all')

                    % Collect all fragments from table
                    selectedTableFragmentIndexes = 1:numel(this.RBEtable.data(1).includedIons);
                    matRad_cfg.dispInfo(sprintf('%d fragments found in RBEtable\n', numel(selectedTableFragmentIndexes)));
                else
                    matRad_cfg.dispError(sprintf('Unrecognized option: "%s" for includedFragments', includedFragments));
                end
            elseif isstruct(includedFragments)
                
                if ~isfield(includedFragments, 'Z')
                    matRad_cfg.dispWarning('No fragment Z selected, setting default');
                    includedFragments.Z = this.defaultFragmentZ;
                end

                if isfield(includedFragments, 'A') && numel(includedFragments.A) == numel(includedFragments.Z) 
                    A = includedFragments.A;
                else
                    A = zeros(size(includedFragments.Z));
                end

                Z = includedFragments.Z;

                selectedTableFragmentIndexes = arrayfun(@(z,a) this.getIndexInTableForFragment(z,a), Z,A,'UniformOutput',false);
                selectedTableFragmentIndexes = cell2mat(selectedTableFragmentIndexes);

            else
                matRad_cfg.dispError('Invalid format for includedFragments');
                selectedTableFragmentIndexes = [];
            end

            this.fragmentIndexesInTable = selectedTableFragmentIndexes;
        end

        % function energies = getFragmentEnergiesFromTable(this, abRatioIndex, fragment, eUnits)
        %     % Converts energies in MeV/u to MeV if necessary. Assumes that
        %     % the base data spectra are given in total energy. TODO: find a
        %     % consisten way of providing this information and check
        %     % consistency between RBEtable and base data.
        % 
        %     energies = this.RBEtable.data(abRatioIndex).energies;
        % 
        %     if strcmp(eUnits, 'MeV/u')
        %         switch fragment
        % 
        %             case 'H1'
        %                 mass = 1;
        % 
        %             case 'C'
        %                 mass = 12;
        % 
        %         end
        % 
        %         energies = energies*mass;
        %     end
        % 
        % end
        % function energies = getFragmentEnergiesFromTable(this,abRatioIndex,fragment, eUnits)
        %     %     % Converts energies in MeV/u to MeV if necessary. Assumes that
        %     %     % the base data spectra are given in total energy. TODO: find a
        %     %     % consisten way of providing this information and check
        %     %     % consistency between RBEtable and base data.
        %     %
        %     energies = this.RBEtable.data(abRatioIndex).energies;
        %     % this.RBEtable.meta.energyDef = input('What unit did you use, MeV/u or MeV?',"s");
        % 
        % 
        %     if strcmp(eUnits, 'MeV/u')
        %         if ~isfield(this.RBEtable.data(1),'convFragEnergy')
        %             v = zeros(1,10);
        %             for i = 1:size(this.RBEtable.data,2)
        %                 this.RBEtable.data(i).convFragEnergy = struct('H1',v,'He',v,'Li',v,'Be',v,'B',v,'C',v);
        %             end
        %         end
        %         %this.RBEtable.meta.consideredFragments = input('Which fragments do you want to look at? Create a string array!');
        %         this.RBEtable.meta.consideredFragments = fragment;
        %         switch this.RBEtable.meta.consideredFragments
        % 
        %             case "H1"
        %                 mass = 1.00794;
        %                 for i = 1: size(this.RBEtable.data,2)
        %                     this.RBEtable.data(i).energies = this.RBEtable.data(i).energies*mass;
        %                     this.RBEtable.data(i).convFragEnergy.H1 = this.RBEtable.data(i).energies;
        %                 end
        %             case "He"
        %                 mass = 4.002602;
        %                 for i = 1: size(this.RBEtable.data,2)
        %                     this.RBEtable.data(i).energies = this.RBEtable.data(i).energies*mass;
        %                     this.RBEtable.data(i).convFragEnergy.He = this.RBEtable.data(i).energies;
        %                 end
        %             case "Li"
        %                 mass = 6.941;
        %                 for i = 1: size(this.RBEtable.data,2)
        %                     this.RBEtable.data(i).energies = this.RBEtable.data(i).energies*mass;
        %                     this.RBEtable.data(i).convFragEnergy.Li = this.RBEtable.data(i).energies;
        %                 end
        %             case "Be"
        %                 mass = 9.012182;
        %                 for i = 1: size(this.RBEtable.data,2)
        %                     this.RBEtable.data(i).energies = this.RBEtable.data(i).energies*mass;
        %                     this.RBEtable.data(i).convFragEnergy.Be = this.RBEtable.data(i).energies;
        %                 end
        %             case "B"
        %                 mass = 10.811;
        %                 for i = 1: size(this.RBEtable.data,2)
        %                     this.RBEtable.data(i).energies = this.RBEtable.data(i).energies*mass;
        %                     this.RBEtable.data(i).convFragEnergy.B = this.RBEtable.data(i).energies;
        %                 end
        %             case "C"
        %                 mass = 12.0107;
        %                 for i = 1: size(this.RBEtable.data,2)
        %                     this.RBEtable.data(i).energies = this.RBEtable.data(i).energies*mass;
        %                     this.RBEtable.data(i).convFragEnergy.C = this.RBEtable.data(i).energies;
        %                 end
        % 
        %         end
        %     end
        % end

        function fragmentIndexesInBaseData = selectBaseDataFragmentsFromMachine(this, machine)
            % This function looks in the base data and checks that the
            % available fluence spectra in the machine are available for
            % the selected fragments

            matRad_cfg = MatRad_Config.instance();

            fragmentIndexesInBaseData = [];

            baseDataFragmentsZ = [machine.data(1).(this.weightBy).spectra.Z];
            baseDataFragmentsA = [machine.data(1).(this.weightBy).spectra.A];

            baseDataZA = [baseDataFragmentsZ', baseDataFragmentsA'];

            tableIncludedFragmentsZ = [this.availableFragmentsInTable(this.fragmentIndexesInTable).Z];
            tableIncludedFragmentsA = [this.availableFragmentsInTable(this.fragmentIndexesInTable).A];

            for fragIdx=1:numel(tableIncludedFragmentsZ)
                currZA = [tableIncludedFragmentsZ(fragIdx), tableIncludedFragmentsA(fragIdx)];

                currBaseDataIndex = find(ismember(baseDataZA, currZA, 'rows'));

                if isempty(currBaseDataIndex)

                    % Check for NaN values in base data
                    if any(currZA(1) == baseDataZA(:,1))
                        matchZ = find(currZA(1) == baseDataZA(:,1));
                        if any(isnan(baseDataZA(matchZ,2)))
                            currBaseDataIndex = find(ismember(baseDataZA(:,1), currZA(1)));
                        end
                    end
                end
                   
                if isempty(currBaseDataIndex) && ~strcmp(machine.meta.radiationMode, 'protons')
                    matRad_cfg.dispError('One or more fragments included in the RBE table are not available in the base data kernel.');
                end
                
                if ~isempty(currBaseDataIndex)
                    fragmentIndexesInBaseData = [fragmentIndexesInBaseData, currBaseDataIndex];
                end
            
            end

            fragmentIndexesInBaseData = unique(fragmentIndexesInBaseData);
            this.fragmentIndexesInBaseData = fragmentIndexesInBaseData;

        end

    end

    methods

        function assignDefaultProperties(this)
            this.defaultRBETable =  'RBEtable_rapidLEMI_testTable';
            %this.tableFragmentIndexes = 1;
            this.defaultFragmentZ = 1;
            this.defaultGenericFragments.Z = [1,2,3,4,5,6,7,8];
            this.defaultGenericFragments.A = [1,4,7,9,11,12,14,16];
        end

        function updateRBEtable(this)
            
            this.RBEtable = this.loadRBEtable(this.RBEtableName);

            this.availableAlphaInTable = [this.RBEtable.data(:).alphaX]';
            this.availableBetaInTable  = [this.RBEtable.data(:).betaX]';
            this.availableFragmentsInTable = [this.RBEtable.data(1).includedIons];

            if ~isempty(this.includedFragments)
                this.updateFragmentIdxInTable(this.includedFragments);
            end
        end

        function tableData = getTableDataForAlphaBeta(this, alphaX, betaX)
            % Scroll the RBE table and check if the requested alpha beta
            % ratio is available

            tableData = [];
            if isempty(this.RBEtable)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('RBEtable not found');
            else
                for i=1:numel(this.RBEtable.data)
                    if (this.RBEtable.data(i).alphaX == alphaX) && (this.RBEtable.data(i).betaX == betaX)
                        tableData = this.RBEtable.data(i);
                        continue;
                    end
                end
                if isempty(tableData)
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispWarning(sprintf('No match found for alpha/beta = %f/%f in RBE table', alphaX, betaX));                    
                end
            end
        end

        function idx = getIndexInTableForFragment(this,Z,A)
            if isempty(this.RBEtable)
                return;
            end

            if ~exist('A', 'var') || isempty(A) || A==0
                A = this.getDefaultAvalueForZ(Z);
            end

            availableZ = [this.RBEtable.data(1).includedIons.Z];
            availableA = [this.RBEtable.data(1).includedIons.A];

            tmpZIdx = find(availableZ == Z);
            tmpAidx = find(availableA == A);

            idx = intersect(tmpZIdx, tmpAidx);

            if numel(idx)>1
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Multiple data for the same ion detected in the RBE table, this should not happen.');
            elseif isempty(idx)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning(sprintf('Unable to find fragment with Z=%d and A=%d in RBE table', Z,A));
            end
        end

        function A = getDefaultAvalueForZ(this, Z)
            idx = find(this.defaultGenericFragments.Z == Z);
            A   = this.defaultGenericFragments.A(idx);
        end

      
    end

    
    methods %(Setters)
        
        function set.RBEtableName(this, value)

            this.RBEtableName = value;
            this.updateRBEtable();

        end

        function set.includedFragments(this, value)

            this.updateFragmentIdxInTable(value);

        end
    
    end

    methods (Static)
        function RBEtable = loadRBEtable(fileName)

            % This function loads the specified RBE table
            matRad_cfg = MatRad_Config.instance();

            searchPath = {fullfile(matRad_cfg.matRadSrcRoot,'bioModels','RBEtables'),...    % default matrad folder
                fullfile(matRad_cfg.primaryUserFolder, 'RBEtables')};             % user defined RBE table

            try
                load(fullfile(searchPath{1}, [fileName, '.mat']), 'RBEtable');

            catch
                try
                    load(fullfile(searchPath{2}, [fileName, '.mat']), 'RBEtable');

                catch
                    matRad_cfg.dispError('Cannot find RBEtable: %s', fileName);
                end
            end
        end



            
        function checkTableConsistency(RBEtable, fragments)

           % This is redundant at the moment
            matRad_cfg = MatRad_Config.instance();

            availableZs = [RBEtable.data(1).includedIons.Z];
            
            if any(~ismember(fragments, availableZs))
                excludedFragments = fragments(~ismember(fragments, availableZs));
                matRad_cfg.dispError('Included RBE table does not contain information %s,',excludedFragments);
            end
          
        end


        function checkRBEtableStructure(RBEtable)
            % Additional function to check structure of the Table

            matRad_cfg = MatRad_Config.instance();

            % Check structure of the table
            if ~isstruct(RBEtable)
                matRad_cfg.dispError('Provided Table is not a struct!');
            end

            if ~isequal(fieldnames(RBEtable), {'meta', 'data'}')

                matRad_cfg.dispError('Provided Table does not contain meta and data fields');

            end

            numOfTissues = numel(RBEtable.data);
            tissueAlphaBetaRatios = [RBEtable.data.alphaX]./[RBEtable.data.betaX];

            includedIons = RBEtable.data(1).includedIons;

            matRad_cfg.dispInfo('%2u tissues found in RBEtable\n', numOfTissues);
            matRad_cfg.dispInfo('with alpha/beta reatios of:\n');

            for tIdx = tissueAlphaBetaRatios
                matRad_cfg.dispInfo('\t %2.3f:\n', tIdx);
            end

            matRad_cfg.dispInfo('Data available for ions: ');
            for fIdx = includedIons
                matRad_cfg.dispInfo('%s ', fIdx{1});
            end

            matRad_cfg.dispInfo('\n');

        end

    end
end