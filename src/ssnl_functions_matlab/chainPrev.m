function chainPrev(this,prev,rules)

rules.fieldsUse = sort(rules.fieldsUse);
[fieldMap,rules] = checkInput(prev,rules);

if rules.props_mixType == 'SHG'
    this.chain_Fields = complex(zeros(2,1,2,size(prev.eField,4)));
elseif rules.props_mixType == 'SFG'
    this.chain_Fields = complex(zeros(2,1,3,size(prev.eField,4)));
end


this.chain_Fields(:,1,fieldMap(:,1),:) =...
    prev.eField(:,end,fieldMap(:,2),:);

this.chain_Fields(:,1,fieldMap(:,1),:) =...
    this.chain_Fields(:,1,fieldMap(:,1),:) .*...
    (prev.props_spotRad/rules.props_spotRad);



if isempty(rules.props_taus)
    
    if length(rules.props_lams) == 2
        prompt = {'Transform Limit of Fund:',...
            'Transform Limit of SH:'};
        defInput = {'330 fs',...
            '20 fs'};
    elseif length(rules.props_lams) == 3
        prompt = {'Transform Limit of Signal:',...
            'Transform Limit of Idler:',...
            'Transform Limit of Pump:'};
        defInput = {'330 fs',...
            '330 fs',...
            '20 fs'};
    end
    inputTitle = 'Specify Pulse Durations';
    dims = [1 100];
    opts.Interpreter = 'tex';
    userInput = inputdlg(prompt,inputTitle,dims,defInput,opts);
    
    % Check if user input something (avoids NaN)
    if ~size(userInput)
        error('Dialog Box Closed. No user input.');
    end
    
    userInput = parseUnits(userInput);
    
    for ii = 1:length(rules.props_lams)
        rules.props_taus(ii) = userInput{ii};
    end
    
end
if isempty(rules.props_energies)
    
    if length(rules.props_lams) == 2
        prompt = {'Energy of Fund:',...
            'Energy of SH:'};
        defInput = {'25 uJ',...
            '0 uJ'};
    elseif length(rules.props_lams) == 3
        prompt = {'Energy of Signal:',...
            'Energy of Idler:',...
            'Energy of Pump:'};
        defInput = {'25 uJ',...
            '25 uJ',...
            '0 uJ'};
    end
    inputTitle = 'Specify Pulse Durations';
    dims = [1 100];
    opts.Interpreter = 'tex';
    userInput = inputdlg(prompt,inputTitle,dims,defInput,opts);
    
    % Check if user input something (avoids NaN)
    if ~size(userInput)
        error('Dialog Box Closed. No user input.');
    end
    
    userInput = parseUnits(userInput);
    
    for ii = 1:length(rules.props_lams)
        rules.props_energies(ii) = userInput{ii};
    end
    
end
if isempty(rules.props_specPhases)
    
    prompt = {'ps or fs?'
        };
    inputTitle = 'Spectral Phase Unit';
    dims = [1 100];
    defInput = {'ps'
        };
    opts.Interpreter = 'tex';
    userInput = inputdlg(prompt,inputTitle,dims,defInput,opts);
    
    % Check if user input something (avoids NaN)
    if ~size(userInput)
        error('Dialog Box Closed. No user input.');
    end
    
    if strcmpi(userInput{1},'ps')
        pFac = 10^-12;
    elseif strcmpi(userInput{1},'ps')
        pFac = 10^-15;
    else
        error('Trying to be sneaky, eh?')
    end
    
    if length(rules.props_lams) == 2
        titles = {'Fund','SH'};
        defGVD = {'0.5','0'};
    elseif length(rules.props_lams) == 3
        titles = {'Signal','Idler','Pump'};
        defGVD = {'0.5','-0.5','0'};
    end
    
    for ii = 1:length(titles)
        if rules.props_energies(ii) ~= 0
            % Get user input to set up beam for GA
            prompt = {'GVD:',...
                'TOD:',...
                '4OD:',...
                '5OD:'
                };
            inputTitle = ['Initial Phase Definition, ',titles{ii}];
            dims = [1 100];
            %         defInput = {'0',...
            defInput = {defGVD{ii},...
                '0',...
                '0',...
                '0'
                };
            opts.Interpreter = 'tex';
            userInput = inputdlg(prompt,inputTitle,dims,defInput,opts);
            
            % Check if user input something (avoids NaN)
            if ~size(userInput)
                error('Dialog Box Closed. No user input.');
            end
            
            userInput = parseUnits(userInput);
            
            rules.props_specPhases(ii,1) = userInput{1}*(pFac^2); %#ok<*ST2NM>
            rules.props_specPhases(ii,2) = userInput{2}*(pFac^3);
            rules.props_specPhases(ii,3) = userInput{3}*(pFac^4);
            rules.props_specPhases(ii,4) = userInput{4}*(pFac^5);
            
        else
            
            rules.props_specPhases(ii,1) = 0;
            rules.props_specPhases(ii,2) = 0;
            rules.props_specPhases(ii,3) = 0;
            rules.props_specPhases(ii,4) = 0;
            
        end
        
    end
    
    
end

if rules.props_taus(1) == 0
    rules.props_taus(1) =...
        this.FWHM(abs(this.chain_Fields(1,1,1,:)).^2) * prev.grid_dt / 10;
end

rules = rmfield(rules,'fieldsUse');
this.inputProperties(rules);

this.props_ks = 2*pi./this.props_lams;
this.props_omegas = this.const_c * this.props_ks;

end

function [fieldMap,rules] = checkInput(prev,rules)

if any(cellfun(@isempty,{rules.fieldsUse,rules.props_mixType rules.props_crys rules.props_len rules.props_theta rules.props_spotRad}))
    error('Missing neccesary properties to chain runs');
end

fieldMap = zeros(length(rules.fieldsUse),2);

if strcmpi(rules.props_mixType,'SFG')
    
%     emptVal = (rules.props_lams == 0 | isnan(rules.props_lams));
    nanVal = isnan(rules.props_lams);
    if ~any(nanVal)
        error('Third wavelength must be specified with NaN');
    end 
    if sum(nanVal) > 1
        error('Not enough wavelengths specified for SFG')
    end
    
    if isequal([0,0,1],nanVal)
        if rules.props_lams(1) <= rules.props_lams(2)
            rules.props_lams(3) = 1./( (1./rules.props_lams(1)) + (1./rules.props_lams(2)) );
            for ii = 1:length(rules.fieldsUse)
                fieldMap(ii,:) = [find(rules.props_lams == prev.props_lams(rules.fieldsUse(ii))) rules.fieldsUse(ii)];
            end
        else
            error('Signal wavelength needs to be less than Idler');
        end
    elseif isequal([0;1;0],nanVal)
        if (rules.props_lams(1) > rules.props_lams(3)) && (rules.props_lams(1) <= 2*rules.props_lams(3))
            rules.props_lams(2) = 1./( (1./rules.props_lams(3)) - (1./rules.props_lams(1)) );
            fieldMap = [[1;3] rules.fieldsUse];
        else
            error('Signal wavelength needs to be greater than and less than twice Pump');
        end
    elseif isequal([1;0;0],nanVal)
        if rules.props_lams(2) >= 2*rules.props_lams(3)
            rules.props_lams(1) = 1./( (1./rules.props_lams(3)) - (1./rules.props_lams(2)) );
            fieldMap = [[2;3] rules.fieldsUse];
        else
            error('Idler wavelength needs to be greater than twice the Pump');
        end
    end
    
    
elseif strcmpi(rules.props_mixType,'SHG')
    
%     emptVal = (rules.props_lams == 0 | isnan(rules.props_lams));
    nanVal = isnan(rules.props_lams);
    if sum(nanVal) > 1
        error('Not enough wavelengths specified for SFG')
    end
    
    if isequal([1;0],nanVal) %#ok<*BDSCA>
        error('Fundamental must be specified for SHG');
    else
        rules.props_lams(2) = rules.props_lams(1)/2;
    end
    
    fieldMap = [1 rules.fieldsUse];
    
    
end


end

function usrOut = parseUnits(usrIn)

usrOut = cell(1,length(usrIn));

for ii = 1:length(usrIn)
    
    matchLst = regexp(usrIn{ii},'-|[0-9]|\.');
    
    numLst = str2num(usrIn{ii}(matchLst)); %#ok<*ST2NM>
    
    if contains(usrIn{ii},'mm') || contains(usrIn{ii},'mJ')
        mFac = 10^-3;
    elseif contains(usrIn{ii},'um') || contains(usrIn{ii},'uJ')
        mFac = 10^-6;
    elseif contains(usrIn{ii},'nm') || contains(usrIn{ii},'nJ') || contains(usrIn{ii},'ns')
        mFac = 10^-9;
    elseif contains(usrIn{ii},'ps')
        mFac = 10^-12;
    elseif contains(usrIn{ii},'fs')
        mFac = 10^-15;
    else
        mFac = 1;
        if length(matchLst) ~= length(usrIn{ii})
            if ~isempty(matchLst)
                warning('Possible bad unit. Assuming it equals 1');
            end
        end
    end
    
    if isempty(matchLst)
        usrOut{ii} = usrIn{ii};
    else
        usrOut{ii} = numLst * mFac;
    end
    
end


end