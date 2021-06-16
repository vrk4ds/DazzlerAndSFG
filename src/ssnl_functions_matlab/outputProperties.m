% Outputs the properties of beamProp
function output = outputProperties(this,str)

% Checks for the properties in this that match 'str'
propList = properties(this);
if ~strcmpi(str,'all')
    propList = propList(contains(propList,str));
end

% Creates a struct to hold the output fields
output = struct;

for ii = 1:length(propList)
    output.(propList{ii}) = this.(propList{ii});
end


end