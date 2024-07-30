function dataOut = matRad_readCsvData(csvFile,cubeDim)
% Read in csv file as table
dataTable = readtable(csvFile,'ReadVariableNames',false);

% this is the number of ReportQuantities contained in that file
numOfReportQuantities = size(dataTable,2)-3;

% Save all scored quantities as cell array and reshape to cubeDim
dataOut = cell(1,numOfReportQuantities);

%Get the indices and do i,j swap
ixi = dataTable.Var2+1;
ixj = dataTable.Var1+1;
ixk = dataTable.Var3+1;
ix = sub2ind(cubeDim,ixi,ixj,ixk);

for i = 1:numOfReportQuantities
    dataOut{i} = zeros(cubeDim);
    dataOut{i}(ix) = dataTable.(['Var' num2str(i+3)]);
end

end