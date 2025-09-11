classdef matRad_TabulatedSpectralKernelBasedModel < matRad_LQRBETabulatedModel
% This is class implementig a spectra-based tabulated RBE model.
% Spectral kernels should be provided in the base data
% The model can handle multiple tissue alphaX/betaX ratio specified by the 
% cst structure, as long as a compatible RBEtable is provided.
%
% Properties of the model that can be defined by the users:
%   RBEtableName        filename of the table to be included
%
%   fragmentIndexesInTable  specifies which fragments to include in the
%                       biological calculation (i.e. 'H', 'He', 'C',...)
%
%   weightBy            specifies which kernel data is used in combination
%                       with the RBE table (i.e. 'Fluence', 'EnergyDeposit')
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
    
    properties (Constant)
        model = 'TAB';
        
        requiredQuantities = {'spectra','physicalDose'};
        kernelQuantities = {'spectra'};
        possibleRadiationModes = {'protons','helium','carbon'};
    end

    properties
        weightBy;
        baseDataKernel;
    end

    methods
        function this = matRad_TabulatedSpectralKernelBasedModel()
            this@matRad_LQRBETabulatedModel();
            this.assignDefaultProperties();
        end


        function bixel = calcBiologicalQuantitiesForBixel(this,bixel,kernel)
        
            % This function implements the variable alpha/beta calculation
            % for the bixel.
            % These values are linearly combined for each fragment and
            % energy. The weighting factor is specified by the selected
            % spectra kernel.
            % Formulation for the mixed field:
            % alpha{d} = sum(sum(alpha(E,F).*Phi(E,F,d))./sum(sum(Phi(E,F,d)));
            % with:
            %   d           radiological depth
            %   E           energy
            %   F           fragment (particle class)
            %   Phi(E,F,d)  Spectum (i.e. Fluence)
            %   alpha(E,F)  alpha value from table for energy E and frgmant F
            % sums are performed over energy and fragments


            bixel = calcBiologicalQuantitiesForBixel@matRad_LQRBETabulatedModel(this,bixel);

            nFragments = numel(this.fragmentIndexesInBaseData);

<<<<<<< HEAD
            % Convert fragments to Z
            for i = 1:numel(this.fragmentsToInclude)
                if strcmp(this.fragmentsToInclude{i}, 'H')
                    includedZ(i) = 1;
                else
                    error('Not implemented')
                end
            end

            % collect the spectra from the bixel, for each fragment
            %spectraEnergies = cellfun(@(fragment) bixel.baseData.Spectra.(this.weightBy).(fragment).energies,this.fragmentsToInclude, 'UniformOutput',false);
            spectraEnergies = arrayfun(@(fragment) bixel.baseData.Fluence(fragment).energyBin,includedZ, 'UniformOutput',false);
            
=======
            % Collect the spectra from the bixel, for each fragment
            % For now assume all fragments have same energy.
            % spectraEnergies = bixel.baseData.spectra.Fluence.energyBin;
             j = 1;
            for i = this.fragmentIndexesInBaseData
                spectraEnergies{1,j} = bixel.baseData.Fluence.spectra(i).energyBin;
                j = j+1;
            end
            % This is to catch situation in which binning in base data
            % starts from zero instead of center of energy bin. This works
            % if bin spacing is linear. If spacing is not linear, provide
            % an energy binning that is bin-centered.
            % if spectraEnergies(1) == 0
            %     spectraEnergies = spectraEnergies + (spectraEnergies(2)-spectraEnergies(1))/2;
            % end
            % 
>>>>>>> dev_quantities_RBE_tabModels_copyRemo
            % Get the tissue classes within the bixel
            bixelTissueIndexes = unique(bixel.vTissueIndex)';
            
            % Interpolate the alpha/beta table for the specific fragment (including primaries)
            % and alphaX/betaX ratio (tissue class)
            % Change here how this is saved into alphaE to have separate
            % fragment fields
            for i=bixelTissueIndexes
                for j = 1:nFragments
                    % [curTissueAlphaE(i,:), curTissueBetaE(i,:)] = cellfun(@(fragment) this.interpolateRBETableForBixel(spectraEnergies, fragment, i),this.fragmentIndexesInTable, 'UniformOutput',false); 
                    % The index of the fragment selected here is taken directly
                    % from the table
                   
                    [curTissueAlphaE, curTissueBetaE] = arrayfun(@(fragment) this.interpolateRBETableForBixel(spectraEnergies, fragment, i),this.fragmentIndexesInTable, 'UniformOutput',false);
                    %dEdx(i,:) = arrayfun(@(fragment) this.interpolateRBETableForBixel(spectraEnergies, fragment,i),this.fragmentIndexesInTable, 'UniformOutput',false);      
                end
            end
             
            % alphaE = arrayfun(@(fragment) horzcat(curTissueAlphaE{:,fragment}), [1:numel(this.fragmentIndexesInTable)], 'UniformOutput',false);
            % betaE  = arrayfun(@(fragment) horzcat(curTissueBetaE{:,fragment}),  [1:numel(this.fragmentIndexesInTable)], 'UniformOutput',false);

            for i = 1:nFragments
                alpha = arrayfun(@(fragment) horzcat(curTissueAlphaE{1,i}), [this.fragmentIndexesInTable(i)], 'UniformOutput',false);
                alphaE{1,i} = alpha{1,1};
                beta = arrayfun(@(fragment) horzcat(curTissueBetaE{1,i}), [this.fragmentIndexesInTable(i)], 'UniformOutput',false);
                betaE{1,i} = beta{1,1};
            end

            % alphaE = cell2struct(alphaE, this.fragmentIndexesInTable, 2);
            % betaE  = cell2struct(betaE,  this.fragmentIndexesInTable, 2);

            % Get the spectra kernels to be used (one for each fragment).
            %bixelSpectra = arrayfun(@(fragmentIdx) squeeze(kernel.(this.weightBy).spectra(fragmentIdx).fluenceSpectrum),this.fragmentIndexesInBaseData, 'UniformOutput',false);
           
            fluenceBixelSpectra = arrayfun(@(fragmentIdx)squeeze(kernel.(this.weightBy).spectra(fragmentIdx).fluenceSpectrum),this.fragmentIndexesInBaseData,'UniformOutput',false);

            % Get normalization for each fragment. This is sum over
            % energies
            % spectraFragmentDenominator = cellfun(@(fragment) sum(bixelSpectra.(fragment),2),this.fragmentIndexesInTable, 'UniformOutput',false);
            spectraFragmentDenominator = cellfun(@(fragmentSpectra) sum(fragmentSpectra,2),fluenceBixelSpectra, 'UniformOutput',false);

            % Get total normalization by summing over the fragments
            spectraBixelDenominator = sum([spectraFragmentDenominator{:}],2);

            % Create containers for alpha and beta for each fragment
            alphaWeightedFragmentSpectra = NaN*ones(numel(bixel.radDepths),nFragments); % only NaN
            betaWeightedFragmentSpectra  = NaN*ones(numel(bixel.radDepths),nFragments);

            for i=bixelTissueIndexes
 
                % Get the spectrum-weighted alpha/beta for each fragment (sums over the energy bins at each bixel radDepth)
                tmpAlphaWeightedFragmentSpectra = arrayfun(@(fragmentIdx) fluenceBixelSpectra{fragmentIdx}*alphaE{fragmentIdx}(:,i), [1:nFragments], 'UniformOutput',false);
                tmpBetaWeightedFragmentSpectra  = arrayfun(@(fragmentIdx) fluenceBixelSpectra{fragmentIdx}*betaE{fragmentIdx}(:,i), [1:nFragments], 'UniformOutput',false);
                %tmpBetaWeightedFragmentSpectra  = cellfun(@(fragment) bixelSpectra.(fragment)*betaE.(fragment)(:,i), this.fragmentIndexesInTable, 'UniformOutput',false);
                
                % Put it in matrix form and sum over the fragments. Only
                % include values for the correct tissue class for each
                % radDepth
                matTmpAlphaWeightedFragmentSpectra = [tmpAlphaWeightedFragmentSpectra{:}];
                matTmpBetaWeightedFragmentSpectra  = [tmpBetaWeightedFragmentSpectra{:}];

                alphaWeightedSpectra(bixel.vTissueIndex == i,:) = sum(matTmpAlphaWeightedFragmentSpectra(bixel.vTissueIndex == i,:),2);
                betaWeightedSpectra(bixel.vTissueIndex == i,:)  = sum(matTmpBetaWeightedFragmentSpectra(bixel.vTissueIndex == i,:),2);
            end

            % Normalize the quantities
            bixel.alpha = alphaWeightedSpectra./spectraBixelDenominator; 
            bixel.beta  = betaWeightedSpectra./spectraBixelDenominator;
        end
    end
    
   methods
        function obj = setFragmentIndexesInBaseData(obj, newVal)
            obj.fragmentIndexesInBaseData = newVal;
        end

        function obj = setTissueAlphaX(obj, newAlphaX)
            obj.tissueAlphaX = newAlphaX;
        end

        function obj = setTissueBetaX(obj, newBetaX)
            obj.tissueBetaX = newBetaX;
        end
    end

    
    methods

        function assignDefaultProperties(this)
            
            assignDefaultProperties@matRad_LQRBETabulatedModel(this);
            this.weightBy = 'Fluence';
            
        end


        function set.weightBy(this, value)

            if ischar(value)
                this.weightBy = value;
                this.updatePropertyValues();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Unrecognized weighting quantity: %s', tostring(value));
            end

        end

       
        function updatePropertyValues(this)
        
            if ~isempty(this.weightBy) && ~isempty(this.fragmentIndexesInTable)
                
                % this.baseDataKernel = cellfun(@(fragment) ['spectra.', this.weightBy, '.', fragment, '.Data'], this.fragmentIndexesInTable, 'UniformOutput',false); %[requiredSpectraData;requriedEnergies];
                % this.baseDataKernel = {this.weightBy}; %cellfun(@(fragment) ['spectra.', this.weightBy], this.fragmentIndexesInTable, 'UniformOutput',false); %[requiredSpectraData;requriedEnergies];
                % this.baseDataKernel = ;
            end


            if ~isempty(this.RBEtable)
                this.checkTableConsistency(this.RBEtable, this.fragmentIndexesInTable);
            end

        end

        function baseDataKernels = getBaseDataKernels(this, baseDataIndexes)
            % baseDataKernels = {sprintf('%s.spectra(%s).fluenceSpectrum', this.weightBy,strjoin(string(baseDataIndexes), ','))};
            baseDataKernels = arrayfun(@(x) sprintf('%s.spectra(%d).fluenceSpectrum',this.weightBy, x), baseDataIndexes, 'UniformOutput', false);
        end

    end
end