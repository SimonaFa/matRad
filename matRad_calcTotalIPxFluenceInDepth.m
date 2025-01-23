function clusterDose = matRad_calcTotalIPxFluenceInDepth(singleData, nameIP, whichZ, whichA)
% This function computes the fluence weighted sum of the IP from the
% fluence spectra of all particles in the machine. Includes all species
% from primaries to secondaries.

%% [1] - Check Fluence Spectra Structure
if ~isfield(singleData, 'Fluence') 
    error('Fluence spectra structure not found.');
else 
    if isfield(singleData.Fluence, 'energyBin')
        energyBin = singleData.Fluence.energyBin;
    else
        error('Energy binning not found.');
    end

    if isfield(singleData.Fluence, 'depthBin')
        depths = singleData.Fluence.depthBin;
    else 
        depths = singleData.depths;
    end
end

%% [2] - Input Parameters and default

typeIP = nameIP(1);
k = str2num(nameIP(2:end));
defZ = [1, 2, 3, 4, 5, 6, 7, 8, -1];
defA = [1, 4, 7, 9, 11, 12, 14, 16, 1];

%% [2] - Allocate Cluster Dose Container

fluenceIP = zeros(size(depths));
if size(fluenceIP, 2) ==1
    fluenceIP = fluenceIP';
end

%% [3] - Loop over Particle Types

if nargin == 2 % take all
    spectraIdx = 1:size(singleData.Fluence.spectra, 2);
elseif nargin == 3
    Zs = [singleData.Fluence.spectra.Z];
    Ms = [singleData.Fluence.spectra.A];
    spectraIdx = [];
    for choice=1:length(whichZ)
        %if whichZ(choice) == 1
        %    spectraIdx(choice) = 1;
        %else
         %   spectraIdx(choice) = find(Zs == whichZ(choice));
        %end
        spectraIdx = [spectraIdx, find(Zs == whichZ(choice))];
    end
elseif nargin ==4
    Zs = [singleData.Fluence.spectra.Z];
    Ms = [singleData.Fluence.spectra.A];
    spectraIdx = [];
    for choice=1:length(whichZ)
        spectraIdx = [spectraIdx, find( (Zs == whichZ(choice)) & (Ms == whichA(choice)) ~= 0) ];
    end
end 

%for particleIdx = 1:size(singleData.fluence.spectra, 2)

for num = 1:length(spectraIdx)
    particleIdx = spectraIdx(num);
    partM = singleData.Fluence.spectra(particleIdx).A;
    if isnan(partM)
        partM = defA(defZ == singleData.Fluence.spectra(particleIdx).Z);
    end
    if singleData.Fluence.spectra(particleIdx).Z == -1 && isfield(singleData.Fluence, 'energyBinEl')
        energyBin = singleData.Fluence.energyBinEl;
    else
    end
    energyPerNucleon = energyBin./partM;
    deltaEPN         = energyPerNucleon(2) - energyPerNucleon(1);
    containerIP      = zeros(size(energyPerNucleon));

    % This will be enlarged for Carbon
    if singleData.Fluence.spectra(particleIdx).Z == 1
        particleName = 'Proton';
    elseif singleData.Fluence.spectra(particleIdx).Z == 2 
        particleName = 'He';
    elseif singleData.Fluence.spectra(particleIdx).Z == 3 
        particleName = 'Li';
    elseif singleData.Fluence.spectra(particleIdx).Z == 4 
        particleName = 'Be';
    elseif singleData.Fluence.spectra(particleIdx).Z == 5 
        particleName = 'B';
    elseif singleData.Fluence.spectra(particleIdx).Z == 6 
        particleName = 'C';
    elseif singleData.Fluence.spectra(particleIdx).Z == 7 
        particleName = 'N';
    elseif singleData.Fluence.spectra(particleIdx).Z == 8
        particleName = 'O';
    elseif singleData.Fluence.spectra(particleIdx).Z == -1
        particleName = 'Electron';
    end
    
    %directoryMCTS           = 'C:\Users\s742o\Work\MCTS-DataBase-UCSF\MCTS-DataBase\MCTS-DataBase';
    fileNameMCTSbaseData    = ['C:\Users\s742o\Work\MCTS-DataBase-UCSF\MCTS-DataBase\MCTS-DataBase', ...
                                '\', particleName, '_', 'IonizationDetail', '_', typeIP, 'k.dat'];

    % This could be substituted with directly reading from the already imported
    % IP table in .mat format.
    IPtableMCTS    = table2array(readtable(fileNameMCTSbaseData));

    % Delicate phase: get the IP vector through interpolation
    for energyBinIdx = 1:length(energyPerNucleon)
        if energyBinIdx == 201
            trythis = 1;
        end
        E1 = energyPerNucleon(energyBinIdx);
        E2 = energyPerNucleon(energyBinIdx) + deltaEPN;

        if E1 >= IPtableMCTS(end, 1)                % I case: Prescription for no intersection, outside boundaries.
            containerIP(energyBinIdx) = IPtableMCTS(end, k+1);
        elseif E2 <= IPtableMCTS(1, 1)
            containerIP(energyBinIdx) = IPtableMCTS(1, k+1);
        else
            % Build the energy and IP comb
            combIdx     = find( IPtableMCTS(:, 1) > E1 & IPtableMCTS(:, 1) < E2 ) ;
            if isempty(combIdx)                     % II case: Inclusion, empty intersection
                containerIP(energyBinIdx) = matRad_interp1(IPtableMCTS(:, 1), IPtableMCTS(:, k+1), (E1+E2)/2);
            else
                energyComb  = [E1; IPtableMCTS(combIdx, 1); E2];
                IPComb      = IPtableMCTS(combIdx, k+1);
                if (E1 <= IPtableMCTS(1, 1)) && (E2 >= IPtableMCTS(end, 1))
                    IP1 = IPtableMCTS(1, k+1);
                    IP2 = IPtableMCTS(end, k+1);
                elseif E1 <= IPtableMCTS(1, 1)          % III case: nonempty intersection, no inclusion
                    IP1 = IPtableMCTS(1, k+1);
                    IP2 = matRad_interp1(IPtableMCTS(:, 1), IPtableMCTS(:, k+1), E2);
                elseif E2 >= IPtableMCTS(end, 1)
                    IP1 = matRad_interp1(IPtableMCTS(:, 1), IPtableMCTS(:, k+1), E1);
                    IP2 = IPtableMCTS(end, k+1);
                else                                % IV case: nonempty intersection, full inclusion 
                    IP1 = matRad_interp1(IPtableMCTS(:, 1), IPtableMCTS(:, k+1), E1);
                    IP2 = matRad_interp1(IPtableMCTS(:, 1), IPtableMCTS(:, k+1), E2);
                end
                IPComb      = [IP1; IPComb; IP2];

                for numEl = 1:length(IPComb)-1
                    containerIP(energyBinIdx) = containerIP(energyBinIdx) + (( IPComb(numEl+1) + IPComb(numEl) )/2)*( ( energyComb(numEl+1) - energyComb(numEl) )/deltaEPN);
                end
            end
        end
    end

    % Contraction
    dataPhi = singleData.Fluence.spectra(particleIdx).fluenceSpectrum;
    %if (singleData.Fluence.spectra(particleIdx).Z == -1)
    %dataPhi(1,:) = zeros(size(singleData.Fluence.spectra(particleIdx).fluenceSpectrum(1,:)));
    %end
    fluenceIP = fluenceIP + containerIP*dataPhi;

end
% Conversion factor from dividing by L = 23.88nm = 23.88 10^-6 (mean chord
% lenght) and rho = 1g/cm^3 = 10^-6 kg/mm^3 (water density).
cf = 10^12/23.88;                % 10^9 becomes 10^7 if the normalization is for cm^2
clusterDose = fluenceIP.*cf;     % CD units [mm^2/kg per primary]

end