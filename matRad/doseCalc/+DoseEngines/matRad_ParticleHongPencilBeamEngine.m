classdef matRad_ParticleHongPencilBeamEngine < DoseEngines.matRad_ParticlePencilBeamEngineAbstract
% matRad_ParticlePencilBeamEngineAbstractGaussian: 
%   Implements an engine for particle based dose calculation 
%   For detailed information see superclass matRad_DoseEngine
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the
% help edit

% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Constant)
           possibleRadiationModes = {'protons', 'helium','carbon'}
           name = 'Hong Particle Pencil-Beam';
           shortName = 'HongPB';

    end
       
    methods 
        
        function this = matRad_ParticleHongPencilBeamEngine(pln)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_ParticleAnalyticalPencilBeamDoseEngine(ct,stf,pln,cst)
            %
            % input
            %   pln:                        matRad plan meta information struct

            if nargin < 1
                pln = [];
            end
             
            this = this@DoseEngines.matRad_ParticlePencilBeamEngineAbstract(pln);
        end
        
    end

    methods (Access = protected)

        function bixel = calcParticleBixel(this,bixel)
            kernels = this.interpolateKernelsInDepth(bixel);
            
            %Lateral Component
            switch this.lateralModel
                case 'single'
                    %compute lateral sigma
                    sigmaSq = kernels.sigma.^2 + bixel.sigmaIniSq;
                    L = exp( -bixel.radialDist_sq ./ (2*sigmaSq))./ (2*pi*sigmaSq);
                case 'double'
                    % compute lateral sigmas
                    sigmaSqNarrow = kernels.sigma1.^2 + bixel.sigmaIniSq;
                    sigmaSqBroad  = kernels.sigma2.^2 + bixel.sigmaIniSq;
    
                    % calculate lateral profile
                    L_Narr =  exp( -bixel.radialDist_sq ./ (2*sigmaSqNarrow))./(2*pi*sigmaSqNarrow);
                    L_Bro  =  exp( -bixel.radialDist_sq./ (2*sigmaSqBroad ))./(2*pi*sigmaSqBroad );
                    L = (1-kernels.weight).*L_Narr + kernels.weight.*L_Bro;
                case 'multi'
                    sigmaSq = kernels.sigmaMulti.^2 + bixel.sigmaIniSq;
                    L = sum([1 - sum(kernels.weightMulti,2), kernels.weightMulti] .* exp(-radialDist_sq ./ (2*sigmaSq))./(2*pi*sigmaSq),2);
                otherwise
                    %Sanity check
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid Lateral Model');
            end
                        
            bixel.physicalDose = bixel.baseData.LatCutOff.CompFac * L .* kernels.Z;
            
            % check if we have valid dose values
            if any(isnan(bixel.physicalDose)) || any(bixel.physicalDose<0)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Error in particle dose calculation.');
            end

            if this.calcLET
                bixel.mLETDose = bixel.physicalDose.*kernels.LET;
            end
            
%{
            if this.calcBioDose
                bixel.mAlphaDose = bixel.physicalDose;
                bixel.mSqrtBetaDose = bixel.physicalDose;
                %From matRad_calcLQParameter
                numOfTissueClass = size(bixel.baseData.alpha,2);
                for i = 1:numOfTissueClass
                    mask = bixel.vTissueIndex == i;
                    if any(mask)
                        bixel.mAlphaDose(mask) = bixel.mAlphaDose(mask) .* kernels.alpha(mask);
                        bixel.mSqrtBetaDose(mask)  = bixel.mSqrtBetaDose(mask) .* kernels.beta(mask);
                    end
                end
            end
%}
            if this.calcBioDose                               
                [bixelAlpha,bixelBeta] = this.bioParam.calcLQParameterForKernel(bixel,kernels);

                bixel.mAlphaDose = bixel.physicalDose .* bixelAlpha;
                bixel.mSqrtBetaDose = bixel.physicalDose .* sqrt(bixelBeta);
            end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if this.calcClusterDose
                if (this.calcClusterDoseFromFluence || ~this.calcCDScatteringFromDose)
                    if isfield(bixel.baseData.Fluence.spectra, 'sigma')
                        error('not implemented \n');
                    elseif isfield(bixel.baseData.Fluence.spectra, 'doubleGauss')
                        error('not implemented \n');
                    elseif isfield(bixel.baseData.Fluence.spectra, 'tripleGauss')

                        % Identify primary particle in order to calculate lateral scattering
                        if strcmp(this.machine.meta.radiationMode, 'carbon')
                            partIdx = 0;
                            for idx = 1:length(bixel.baseData.Fluence.spectra)
                                if bixel.baseData.Fluence.spectra(idx).Z == 6
                                    partIdx = idx;
                                end
                            end
                        elseif strcmp(this.machine.meta.radiationMode, 'helium')
                            partIdx = 0;
                            for idx = 1:length(bixel.baseData.Fluence.spectra)
                                if bixel.baseData.Fluence.spectra(idx).Z == 2
                                    partIdx = idx;
                                end
                            end
                        elseif strcmp(this.machine.meta.radiationMode, 'protons')
                            partIdx = 0;
                            for idx = 1:length(bixel.baseData.Fluence.spectra)
                                if (bixel.baseData.Fluence.spectra(idx).Z == 1) && (bixel.baseData.Fluence.spectra(idx).A == 1)
                                    partIdx = idx;
                                end
                            end
                        else
                            error('primary particle not found \n');
                        end

                        %{
                    sigmaPrimary(partIdx).s1 = sqrt(kernels.Fluence(partIdx).sigma1.^2 + bixel.sigmaIniSq);
                    sigmaPrimary(partIdx).s2 = sqrt(kernels.Fluence(partIdx).sigma2.^2 + bixel.sigmaIniSq);
                    sigmaPrimary(partIdx).s3 = sqrt(kernels.Fluence(partIdx).sigma3.^2 + bixel.sigmaIniSq);

                    LcDPrim = 0;
                    LcDTotal = 0;
                    LcDSec = 0;
                    norm = 0;
                        %}
                        %{
                    for partIdx = 1:numel(kernels.fluence)
                        Lcd = Lcd + kernels.fluence(partIdx).cumFluence .* gaussDist3(sqrt(bixel.radialDist_sq), sigmaFl(partIdx).s1, sigmaFl(partIdx).s2, sigmaFl(partIdx).s3, kernels.fluence(1).w2, kernels.fluence(1).w3);
                        norm = norm + kernels.fluence(partIdx).cumFluence;
                    end
                        %}
                        %Lcd = Lcd ./ (norm);

                        %LcDPrim = gaussDist3(sqrt(bixel.radialDist_sq), sigmaPrimary(partIdx).s1, sigmaPrimary(partIdx).s2, sigmaPrimary(partIdx).s3, kernels.Fluence(partIdx).w2, kernels.Fluence(partIdx).w3);

                        % Calc clusterDose Primary


                        % compute lateral sigmas
                        %sigma1cd = kernels.sigma1.^2 + bixel.sigmaIniSq;
                        %sigma2cd = kernels.sigma2.^2 + bixel.sigmaIniSq;

                        % calculate lateral profile
                        %L_Narr =  exp( -bixel.radialDist_sq ./ (2*sigmaSqNarrow))./(2*pi*sigmaSqNarrow);
                        %L_Bro  =  exp( -bixel.radialDist_sq./ (2*sigmaSqBroad ))./(2*pi*sigmaSqBroad );
                        %L = (1-kernels.weight).*L_Narr + kernels.weight.*L_Bro;
                    end
                end

                % Calc cluster dose Primary
                if this.calcClusterDoseFromFluence
                    clusterDoseTotal = zeros(size(kernels.clusterDoseParticles(1).clusterDoseProfile));
                    if this.calcPrimary
                        clusterDosePrimary = zeros(size(kernels.clusterDoseParticles(1).clusterDoseProfile));
                    end
                    if this.calcSecondary
                        clusterDoseSecondary = zeros(size(kernels.clusterDoseParticles(1).clusterDoseProfile));
                    end
                    for pIdx = 1:(length(bixel.baseData.Fluence.spectra)-1)
                        % Calc sigmas
                        sigmaParticle(pIdx).s1 = sqrt(kernels.Fluence(pIdx).sigma1.^2 + bixel.sigmaIniSq);
                        sigmaParticle(pIdx).s2 = sqrt(kernels.Fluence(pIdx).sigma2.^2 + bixel.sigmaIniSq);
                        sigmaParticle(pIdx).s3 = sqrt(kernels.Fluence(pIdx).sigma3.^2 + bixel.sigmaIniSq);

                        % Lateral scattering
                        LcD = this.tripleGaussRadial( sqrt(bixel.radialDist_sq), sigmaParticle(pIdx).s1, sigmaParticle(pIdx).s2, sigmaParticle(pIdx).s3, kernels.Fluence(partIdx).w2, kernels.Fluence(partIdx).w3 );

                        % calc Cluster Dose
                        clusterDoseTotal = clusterDoseTotal + LcD .* kernels.clusterDoseParticles(pIdx).clusterDoseProfile;
                        if (this.calcPrimary || this.calcSecondary)
                            if (pIdx == partIdx) % Primary
                                clusterDosePrimary = LcD .* kernels.clusterDoseParticles(pIdx).clusterDoseProfile;
                            else
                                clusterDoseSecondary = clusterDoseSecondary + LcD .* kernels.clusterDoseParticles(pIdx).clusterDoseProfile;
                            end
                        end
                    end

                    % Calc Total Cluster Dose
                    bixel.mClusterDose          = clusterDoseTotal;
                    if this.calcPrimary
                        bixel.mClusterDosePrimary   = clusterDosePrimary;
                    end
                    if this.calcSecondary
                        bixel.mClusterDoseSecondary = clusterDoseSecondary;
                    end
                else

                    if isfield(kernels, 'cDoseSigma1') && ~this.calcCDScatteringFromDose
                        % Calc sigmas
                        sigmaParticle.s1 = sqrt(kernels.cDoseSigma1.^2 + bixel.sigmaIniSq);
                        sigmaParticle.s2 = sqrt(kernels.cDoseSigma2.^2 + bixel.sigmaIniSq);
                        sigmaParticle.s3 = sqrt(kernels.cDoseSigma3.^2 + bixel.sigmaIniSq);
                        % Lateral triple Gaussian from Primary
                        LcD = this.tripleGaussRadial( sqrt(bixel.radialDist_sq), sigmaParticle.s1, sigmaParticle.s2, sigmaParticle.s3, kernels.cDoseWeight2, kernels.cDoseWeight3 );
                        bixel.mClusterDose  = LcD .* kernels.clusterDose;
                    else
                        bixel.mClusterDose  = L .* kernels.clusterDose;
                    end

                end
                %{
                bixel.mClusterDose = LcDPrim .* kernels.clusterDose;
                if this.calcPrimary
                    bixel.mClusterDosePrimary   = LcDPrim .* kernels.clusterDosePrimary;
                end
                if this.calcSecondary
                    bixel.mClusterDoseSecondary = LcDPrim .* kernels.clusterDoseSecondary;
                end
                %}
            end

        end
    
        function dist3 = tripleGaussRadial(~, radialDist, sigma1, sigma2, sigma3, weight2, weight3)
            gaussDist = @(r, s) exp( -r.^2 ./ (2*s.^2)) ./ (2*pi*s.^2);
            gaussDist3 = @(r, s1, s2, s3, w2, w3) (1 - w2 - w3) .* gaussDist(r, s1) + w2 .* gaussDist(r, s2) + w3 .* gaussDist(r, s3);
            dist3 = gaussDist3(radialDist, sigma1, sigma2, sigma3, weight2, weight3);
        end
    end

    methods (Static)
        function [available,msg] = isAvailable(pln,machine)   
            % see superclass for information
            
            msg = [];
            available = false;

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            %checkBasic
            try
                checkBasic = isfield(machine,'meta') && isfield(machine,'data');

                %check modality
                checkModality = any(strcmp(DoseEngines.matRad_ParticleHongPencilBeamEngine.possibleRadiationModes, machine.meta.radiationMode));
                
                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end

            checkMeta = all(isfield(machine.meta,{'SAD','BAMStoIsoDist','LUT_bxWidthminFWHM','dataType'}));

            dataType = machine.meta.dataType;
            if strcmp(dataType,'singleGauss')
                checkData = all(isfield(machine.data,{'energy','depths','Z','peakPos','sigma','offset','initFocus'}));
            elseif strcmp(dataType,'doubleGauss')
                checkData = all(isfield(machine.data,{'energy','depths','Z','peakPos','weight','sigma1','sigma2','offset','initFocus'}));
            elseif strcmp(dataType,'multipleGauss')
                checkData = all(isfield(machine.data,{'energy','depths','Z','peakPos','weightMulti','sigmaMulti','offset','initFocus'}));
            else
                checkData = false;
            end
            
            available = checkMeta && checkData;
        end
    end
end

