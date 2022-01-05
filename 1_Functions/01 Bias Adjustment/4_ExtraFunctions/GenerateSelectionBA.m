function selection = GenerateSelectionBA(timefiles)
%GENERATESELECTIONBA Summary of this function goes here
%   Detailed explanation goes here

%% Set-up
filenumber = length(timefiles);
selection = cell(filenumber, 4);

%% Generation

for files=1:filenumber
    selection{files,1} = find(timefiles{files}(:,2) == 12 | timefiles{files}(:,2) == 1 | timefiles{files}(:,2) == 2);
    selection{files,2} = find(timefiles{files}(:,2) == 3 | timefiles{files}(:,2) == 4 | timefiles{files}(:,2) == 5);
    selection{files,3} = find(timefiles{files}(:,2) == 6 | timefiles{files}(:,2) == 7 | timefiles{files}(:,2) == 8);
    selection{files,4} = find(timefiles{files}(:,2) == 9 | timefiles{files}(:,2) == 10 | timefiles{files}(:,2) == 11);
end

