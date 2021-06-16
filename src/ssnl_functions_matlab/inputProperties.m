% Inputs the properties defined by the structure vec
function inputProperties(this,vec)

% Grabs the field names from the vec struct
fieldList = fieldnames(vec);
fieldList = fieldList(~contains(fieldList,'const'));

% This will throw an error if there is not a property of this that is set
for ii = 1:length(fieldList)
    this.(fieldList{ii}) = vec.(fieldList{ii});
end


end