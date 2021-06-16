function setupGUI(this)

% Get type of mixing (3 or 2 wave)
prompt = {'SFG or SHG?'};
inputTitle = 'Choose type of mixing';
dims = [1 100];
defInput = {'SHG'};
opts.Interpreter = 'tex';
userInput = inputdlg(prompt,inputTitle,dims,defInput,opts);

% Check if user input something (avoids NaN)
if ~size(userInput)
    error('Dialog Box Closed. No user input.');
end


if strcmpi(userInput{1},'SFG')
    this.props_mixType = 'SFG';
    this.props_lams = zeros(3,1);
    this.props_omegas = zeros(3,1);
    this.props_ks = zeros(3,1);
    this.props_taus = zeros(3,1);
    this.props_energies = zeros(3,1);
    this.props_specPhases = zeros(3,4);
    mixType = 1;
else
    this.props_mixType = 'SHG';
    this.props_lams = zeros(2,1);
    this.props_omegas = zeros(2,1);
    this.props_ks = zeros(2,1);
    this.props_taus = zeros(2,1);
    this.props_energies = zeros(2,1);
    this.props_specPhases = zeros(2,4);
    mixType = 2;
end

if mixType == 1
    opts.Interpreter = 'tex';
    opts.Default = 'Sig & Idler';
    userInput = questdlg('Which wavelengths will you specify? (SPACE = Return)',...
        'Choose Mixing Type',...
        'Sig & Idler',...
        'Sig & Pump',...
        'Idler & Pump',...
        opts);
    
    switch userInput
        case 'Sig & Idler'
            sfgType = {'Signal','Idler','1030','1030',1};
        case 'Sig & Pump'
            sfgType = {'Signal','Pump','1030','515',2};
        case 'Idler & Pump'
            sfgType = {'Idler','Pump','1030','515',3};
        otherwise
            sfgType = {'Signal','Idler','1030','515',1};
    end
    
end


% Get user input to set up beam for GA

if mixType == 1
    prompt = {['Wavelength of ',sfgType{1},':'],...
        ['Wavelength of ',sfgType{2},':'],...
        'Transform Limit of Signal:',...
        'Transform Limit of Idler:',...
        'Transform Limit of Pump:',...
        'Energy of Signal:',...
        'Energy of Idler:',...
        'Energy of Pump:',...
        'Radius of Spot:'
        };
    defInput = {[sfgType{3},' nm'],...
        [sfgType{4},' nm'],...
        '330 fs',...
        '330 fs',...
        '20 fs',...
        '25 uJ',...
        '25 uJ',...
        '0 uJ',...
        '400 um'
        };
else
    prompt = {'Wavelength of Fund:',...
        'Transform Limit of Fund:',...
        'Transform Limit of SH',...
        'Energy of Fund:',...
        'Energy of SH:',...
        'Radius of Spot:'
        };
    defInput = {'1030 nm',...
        '330 fs',...
        '20 fs',...
        '25 uJ',...
        '0 uJ',...
        '400 um'
        };
end

inputTitle = 'Initial Beam Definition, Include units (ie. ''330 fs'')';
dims = [1 100];
opts.Interpreter = 'tex';
userInput = inputdlg(prompt,inputTitle,dims,defInput,opts);

% Check if user input something (avoids NaN)
if ~size(userInput)
    error('Dialog Box Closed. No user input.');
end

if mixType == 1
    
    % If the tau on field 3 is zero you get 1/0 errors
    if userInput{5} == 0
        userInput{5} = defInput{5};
    end
    
    userInput = parseUnits(userInput);
    
    if sfgType{5} == 1 && userInput{1} > userInput{2}
        error('Signal wavelength must be smaller than Idler');
    elseif sfgType{5} == 2 && ((userInput{1} < userInput{2}) || (userInput{1} > 2*userInput{2}))
        error('Signal wavelength must be larger than and smaller than twice PumpPump');
    elseif sfgType{5} == 3 && userInput{1} <= 2*userInput{2}
        error('Idler wavelength must be larger than twice the Pump');
    end
    
    if sfgType{5} == 1
    this.props_lams(1) = userInput{1}; %#ok<*ST2NM>
    this.props_lams(2) = userInput{2};
    this.props_lams(3) = 1/...
        ( (1/this.props_lams(1)) + (1/this.props_lams(2)) );
    elseif sfgType{5} == 2
    this.props_lams(1) = userInput{1};
    this.props_lams(3) = userInput{2};
    this.props_lams(2) = 1/...
        ( (1/this.props_lams(3)) - (1/this.props_lams(1)) );
    elseif sfgType{5} == 3
    this.props_lams(2) = userInput{1};
    this.props_lams(3) = userInput{2};
    this.props_lams(1) = 1/...
        ( (1/this.props_lams(3)) - (1/this.props_lams(2)) );
    end
    
    this.props_taus(1) = userInput{3};
    this.props_taus(2) = userInput{4};
    this.props_taus(3) = userInput{5};
    this.props_energies(1) = userInput{6};
    this.props_energies(2) = userInput{7};
    this.props_energies(3) = userInput{8};
    this.props_spotRad = userInput{9};
    
    
    this.props_ks = 2*pi./this.props_lams;
    this.props_omegas = this.const_c * this.props_ks;
    
else
    
    % If the tau on field 2 is zero you get 1/0 errors
    if userInput{3} == 0
        userInput{3} = defInput{3};
    end
    
    userInput = parseUnits(userInput);
    
    this.props_lams(1) = userInput{1}; %#ok<*ST2NM>
    this.props_lams(2) = this.props_lams(1) / 2;
    this.props_taus(1) = userInput{2};
    this.props_taus(2) = userInput{3};
    this.props_energies(1) = userInput{4};
    this.props_energies(2) = userInput{5};
    this.props_spotRad = userInput{6};
    
    
    this.props_ks = 2*pi./this.props_lams;
    this.props_omegas = this.const_c * this.props_ks;
    
end


% Get user input to set up beam for GA
prompt = {'Type of Crystal:',...
    'Length of Crystal:',...
    'Angle of Crystal (degrees):'
    };
defInput = {'BBO',...
    '1.75 mm',...
    '23.29'
    };
inputTitle = 'Initial Crystal Definition, Include units (ie. ''330 fs'')';
dims = [1 100];
opts.Interpreter = 'tex';
userInput = inputdlg(prompt,inputTitle,dims,defInput,opts);

% Check if user input something (avoids NaN)
if ~size(userInput)
    error('Dialog Box Closed. No user input.');
end

userInput = parseUnits(userInput);

this.props_crys = userInput{1}; %#ok<*ST2NM>
this.props_len = userInput{2};
this.props_theta = userInput{3};


% Get user input to set up beam for GA
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

if mixType == 1
    titles = {'Signal','Idler','Pump'};
    
    defGVD = {'0.5','-0.5','0'};
    
else
    titles = {'Fund','SH'};
    
    defGVD = {'0.5','0'};
    
end


for ii = 1:length(titles)
    if this.props_energies(ii) ~= 0
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
        
        this.props_specPhases(ii,1) = userInput{1}*(pFac^2); %#ok<*ST2NM>
        this.props_specPhases(ii,2) = userInput{2}*(pFac^3);
        this.props_specPhases(ii,3) = userInput{3}*(pFac^4);
        this.props_specPhases(ii,4) = userInput{4}*(pFac^5);
        
    else
        
        this.props_specPhases(ii,1) = 0;
        this.props_specPhases(ii,2) = 0;
        this.props_specPhases(ii,3) = 0;
        this.props_specPhases(ii,4) = 0;
        
    end
    
end


end

function usrOut = parseUnits(usrIn)

usrOut = cell(1,length(usrIn));

for ii = 1:length(usrIn)
    
    matchLst = regexp(usrIn{ii},'-|[0-9]|\.');
    
    numLst = str2num(usrIn{ii}(matchLst));
    
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