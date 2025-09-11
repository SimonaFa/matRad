classdef matRad_StoppingPowerModel < matRad_LQRBETabulatedModel
% This is class implementig a stopping power to the tabulated RBE model.
% dE/dx data in the RBE table should be provided
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
% Copyright 2025 the matRad development team.
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
        model = 'StPo';
        
        requiredQuantities = {'spectra','physicalDose'};
        kernelQuantities = {'spectra'};
        possibleRadiationModes = {'protons','helium','carbon'};
    end

    properties
        weightBy;
        baseDataKernel;
    end

    methods
        function this = matRad_StoppingPowerModel()
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

            % Collect the spectra from the bixel, for each fragment
            % For now assume all fragments have same energy.
            % spectraEnergies = bixel.baseData.spectra.Fluence.energyBin;
            j = 1;
            for i = this.fragmentIndexesInBaseData
                spectraEnergies{1,j} = bixel.baseData.Fluence.spectra(i).energyBin;
                j = j+1;
            end
            
            % Get the tissue classes within the bixel
            bixelTissueIndexes = unique(bixel.vTissueIndex)';
            
            % Interpolate the alpha/beta table for the specific fragment (including primaries)
            % and alphaX/betaX ratio (tissue class)
            % Change here how this is saved into alphaE to have separate
            % fragment fields
            for i=bixelTissueIndexes
                %for j = 1:nFragments
                    % [curTissueAlphaE(i,:), curTissueBetaE(i,:)] = cellfun(@(fragment) this.interpolateRBETableForBixel(spectraEnergies, fragment, i),this.fragmentIndexesInTable, 'UniformOutput',false); 
                    % The index of the fragment selected here is taken directly
                    % from the table
                   
                    [curTissueAlphaE{1,i}, curTissueBetaE{1,i}] = arrayfun(@(fragment) this.interpolateRBETableForBixel(spectraEnergies, fragment, i),this.fragmentIndexesInTable, 'UniformOutput',false);
                    %dEdx(i,:) = arrayfun(@(fragment) this.interpolateRBETableForBixel(spectraEnergies, fragment,i),this.fragmentIndexesInTable, 'UniformOutput',false);      
                %end
                % alpha{1,i} = arrayfun(@(fragment) horzcat(curTissueAlphaE{i}{fragment}{:}), this.fragmentIndexesInTable, 'UniformOutput', false);
                alphaE{1,i} = curTissueAlphaE{i};
                % %alphaE{1,i} = alpha{1,1};
                % beta{1,i} = arrayfun(@(fragment) horzcat(curTissueBetaE{i}{fragment}{:}), [this.fragmentIndexesInTable], 'UniformOutput',false);
                %betaE{1,i} = beta{1,1};
                betaE{1,i} = curTissueBetaE{i};
            end


            for i = 1:nFragments                
                dEdx{1,i} = matRad_interp1(this.RBEtable.data(bixelTissueIndexes).energies, this.RBEtable.data(bixelTissueIndexes).dEdx(:,i), spectraEnergies{1,i});
                % dEdx = arrayfun(@(fragment) horzcat(dEdx{1,i}), [this.fragmentIndexesInTable(i)], 'UniformOutput',false);
                % dEdxE{1,i} = dEdx{1,1};
            end

            % for i = 1:numel(this.fragmentIndexesInBaseData)
            %     fragIdx = this.fragmentIndexesInBaseData(i);
            %     meanE{i} = spectraEnergies .* bixel.baseData.Fluence.spectra(fragIdx).fluenceSpectrum(:,2:end);
            %     dE{i} = meanE{i}(2:end) - meanE{i}(1:end-1);
            %     dx{i} = (bixel.baseData.depths(3:end) - bixel.baseData.depths(2:end-1))';
            %     dEdX{i} = interp1(meanE{i}(1:end-1), -dE{i}./dx{i},spectraEnergies,'spline','extrap');
            % end
            
            % Get the spectra kernels to be used (one for each fragment).
            %bixelSpectra = arrayfun(@(fragmentIdx) squeeze(kernel.(this.weightBy).spectra(fragmentIdx).fluenceSpectrum),this.fragmentIndexesInBaseData, 'UniformOutput',false);
           
            fluenceBixelSpectra = arrayfun(@(fragmentIdx)squeeze(kernel.(this.weightBy).spectra(fragmentIdx).fluenceSpectrum),this.fragmentIndexesInBaseData,'UniformOutput',false);
            %bixelSpectra = arrayfun(@(fragmentIdx)fluenceBixelSpectra{fragmentIdx}.*dEdX{fragmentIdx},1:numel(fluenceBixelSpectra),'UniformOutput',false);
            bixelSpectra = arrayfun(@(fragmentIdx)fluenceBixelSpectra{fragmentIdx}.*dEdx{fragmentIdx}',1:numel(fluenceBixelSpectra),'UniformOutput',false); % fluence times stopping power to get dose + the sum
            
            % Get normalization for each fragment. This is sum over
            % energies
            % spectraFragmentDenominator = cellfun(@(fragment) sum(bixelSpectra.(fragment),2),this.fragmentIndexesInTable, 'UniformOutput',false);
            %spectraFragmentDenominator = cellfun(@(fragmentSpectra) sum(fragmentSpectra,2),bixelSpectra, 'UniformOutput',false);
            
            spectraFragmentDenominator = sum(cell2mat(bixelSpectra), 2); % sum over the fragments

            % Get total normalization by summing over the fragments
            %spectraBixelDenominator = sum([spectraFragmentDenominator{:}],2);

            % Create containers for alpha and beta for each fragment
            alphaWeightedFragmentSpectra = NaN*ones(numel(bixel.radDepths),nFragments); % only NaN
            betaWeightedFragmentSpectra  = NaN*ones(numel(bixel.radDepths),nFragments);

            for i=bixelTissueIndexes
 
                % Get the spectrum-weighted alpha/beta for each fragment (sums over the energy bins at each bixel radDepth)
                tmpAlphaWeightedFragmentSpectra = arrayfun(@(fragmentIdx) bixelSpectra{fragmentIdx}*(alphaE{1,bixelTissueIndexes}{1,fragmentIdx}), [1:nFragments], 'UniformOutput',false);
                tmpBetaWeightedFragmentSpectra  = arrayfun(@(fragmentIdx) bixelSpectra{fragmentIdx}*(betaE{1,bixelTissueIndexes}{1,fragmentIdx}), [1:nFragments], 'UniformOutput',false);
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
            bixel.alpha = alphaWeightedSpectra./spectraFragmentDenominator; 
            bixel.beta  = betaWeightedSpectra./spectraFragmentDenominator;
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