% matRad script
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_rc

% load patient data, i.e. ct, voi, cst

%load HEAD_AND_NECK
%load TG119.mat
%load PROSTATE_simultaneous5_cdDom_doseConstr.mat
%load PROSTATE_dominantCD.mat
%load PROSTATE_new2.mat %GOOD
%load PROSTATE_new1_sim5.mat
%load LIVER_new_sim5.mat % last
%load BOXPHANTOM_simultaneous4.mat
%load BOXPHANTOM_mod_onlytarget.mat
%load BOXPHANTOM_simultaneous4.mat
%load BOXPHANTOM_new_sim1.mat
%load PROSTATE.mat
%load BOXPHANTOM_dominantCD.mat
%load PROSTATE_new4_sim1.mat
%load PROSTATE_new8_VIII.mat
%load BOXPHANTOM.mat

% meta information for treatment plan

pln.radiationMode   = 'protons';     % either photons / protons / carbon%
%pln.radiationMode    = 'helium';
%pln.machine         = 'Generic';
%pln.machine         = 'generic_TOPAS_clusterdose_draft';
%pln.machine         = 'generic_cd_draft';
%pln.machine         = 'generic_TOPAS_500_cd';
%pln.machine         = 'generic_TOPAS_cd';
%pln.machine         = 'Generic_fluence_cluster_Facchiano';

pln.machine         = 'Generic_fluence_cluster_Facchiano_interp.mat';
%pln.machine         = 'Generic_Fluence_Cluster_Facchiano_15.11.2023_focusmod';

pln.numOfFractions  = 30;

% beam geometry settings
pln.propStf.bixelWidth      = 5; %1000 [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [-90]; % [?]
pln.propStf.couchAngles     = [0]; % [?]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propStf.longitudinalSpotSpacing = 5; %20;

% dose calculation settings
%pln.propDoseCalc.engine = 'none';
%pln.propDoseCalc.engine  = 'Analytical Bragg Peak';
pln.propDoseCalc.engine  = 'Particle Pencil-Beam';

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]
%pln.propDoseCalc.doseGrid.resolution = ct.resolution;

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none';
                                      % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles

pln.propSeq.runSequencing   = true;  % true: run sequencing, false: don't / will be ignored for particles and also triggered by runDAO below

%% initial visualization and change objective function settings if desired
matRadGUI

%% generate steering file
%stf = matRad_generateSingleBixelStf(ct,cst,pln);
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
dij = matRad_calcDose(ct,cst,stf,pln);

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
%{
resultGUI.physicalDose = reshape(dij.physicalDose{1}*ones(dij.totalNumOfBixels,1), dij.doseGrid.dimensions);
for i = 1:dij.totalNumOfBixels
    weight = zeros(dij.totalNumOfBixels,1);
    weight(i) = 1;
    resultGUI.(['physicalDose_' num2str(i)]) = reshape(dij.physicalDose{1}*weight, dij.doseGrid.dimensions);
end
%}
%% sequencing
resultGUI = matRad_sequencing(resultGUI,stf,dij,pln);


%% DAO
if strcmp(pln.radiationMode,'photons') && pln.propOpt.runDAO
   resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
   matRad_visApertureInfo(resultGUI.apertureInfo);
end

%% start gui for visualization of result
matRadGUI

%% indicator calculation and show DVH and QI
[dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI);

