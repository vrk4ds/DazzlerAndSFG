function plotField(this, fInd, fNum, range, norm)

if ~exist('norm','var') || isempty(norm)
    norm = 0;
end


checkInputs(this,fInd,fNum,range,norm)

pInd = cell(1,4);

pInd(1:3) = parsePlotType(fInd,fNum,this.grid_nZSteps);
pInd(4) = findRange(this,pInd,range);


if norm == 0
    
    if length(pInd{2}) == 1 %% 2D %%
        
        if pInd{1} == 1
            
            plot(...
                this.list_t(pInd{4})/10^-12,...
                abs(squeeze( this.eField(pInd{1,:}) )).^2 ...
                );
            xlim([...
                this.list_t(pInd{4}(1))/10^-12 ...
                this.list_t(pInd{4}(end))/10^-12 ...
                ])
            
            xlabel('Time (ps)','Interpreter','latex')
            ylabel('Intesity ($\frac{W}{m^2}$)','Interpreter','latex')
            
        elseif pInd{1} == 2
            
            plot(...
                this.list_lambda(pInd{3},pInd{4})/10^-9,...
                abs(squeeze( this.eField(pInd{1,:}) )).^2 ...
                );
            xlim([...
                this.list_lambda(pInd{3},pInd{4}(end))/10^-9 ...
                this.list_lambda(pInd{3},pInd{4}(1))/10^-9 ...
                ])
            
            xlabel('Wavelength (nm)','Interpreter','latex')
            ylabel('Intesity ($\frac{W}{m^2 Hz}$)','Interpreter','latex')
            
        end
        
        
    elseif length(pInd{2}) > 1 %% 3D %%
        
        if pInd{1} == 1
            
            xVals = (pInd{2}/max(pInd{2}))*this.props_len/10^-3;
            tVals = this.list_t(pInd{4})/10^-12;
            
            surf( xVals, tVals,...
                abs(squeeze( this.eField(pInd{:}) )').^2,...
                'EdgeColor','none','FaceLighting','gouraud');
%             view([45 30])
             view([0 90])
            
            xlim([xVals(1) xVals(end)])
            ylim([tVals(1) tVals(end)])
            xticks(unique(round(linspace(xVals(1),xVals(end),max(ceil(length(xVals)/8),5)),2)))
            
            xlabel('Z-Coordinate (mm)','Interpreter','latex')
            ylabel('Time (ps)','Interpreter','latex')
            zlabel('Intesity ($\frac{W}{m^2}$)','Interpreter','latex')
            
        elseif pInd{1} == 2
            
            xVals = (pInd{2}/max(pInd{2}))*this.props_len/10^-3;
            tVals = this.list_lambda(pInd{3},pInd{4})/10^-9;
            
            surf( xVals, tVals,...
                abs(squeeze( this.eField(pInd{:}) )').^2,...
                'EdgeColor','none','FaceLighting','gouraud');
%             view([45 30])
             view([0 90])
            
            xlim([xVals(1) xVals(end)])
            ylim([tVals(end) tVals(1)])
            xticks(unique(round(linspace(xVals(1),xVals(end),max(ceil(length(xVals)/8),5)),2)))
            
            xlabel('Z-Coordinate (mm)','Interpreter','latex')
            ylabel('Wavelength (nm)','Interpreter','latex')
            zlabel('Intesity ($\frac{W}{m^2 Hz}$)','Interpreter','latex')
            
        end
                
    end
    
else
    
    if length(pInd{2}) == 1 %% 2D %%
        
        if pInd{1} == 1
            
            plot(...
                this.list_t(pInd{4})/10^-12,...
                abs(squeeze( this.eField(pInd{1,:}) )).^2 /...
                max(abs(squeeze( this.eField(pInd{1,:}) )).^2)...
                );
            xlim([...
                this.list_t(pInd{4}(1))/10^-12 ...
                this.list_t(pInd{4}(end))/10^-12 ...
                ])
            ylim([0 1])
            
            xlabel('Time (ps)','Interpreter','latex')
            
        elseif pInd{1} == 2
            
            plot(...
                this.list_lambda(pInd{3},pInd{4})/10^-9,...
                abs(squeeze( this.eField(pInd{1,:}) )).^2 /...
                max(abs(squeeze( this.eField(pInd{1,:}) )).^2)...
                );
            
            xlim([...
                this.list_lambda(pInd{3},pInd{4}(end))/10^-9 ...
                this.list_lambda(pInd{3},pInd{4}(1))/10^-9 ...
                ])
            ylim([0 1])
            
            xlabel('Wavelength (nm)','Interpreter','latex')
            
        end
        
        ylabel('Intesity (Norm.)','Interpreter','latex')
        
    elseif length(pInd{2}) > 1 %% 3D %%
        
        if pInd{1} == 1
            
            xVals = (pInd{2}/max(pInd{2}))*this.props_len/10^-3;
            tVals = this.list_t(pInd{4})/10^-12;
            
            surf( xVals, tVals,...
                abs(squeeze( this.eField(pInd{:}) )').^2 /...
                max(max(abs(squeeze( this.eField(pInd{:}) )').^2)),...
                'EdgeColor','none','FaceLighting','gouraud');
%             view([45 30])
             view([0 90])
            
            xlim([xVals(1) xVals(end)])
            ylim([tVals(1) tVals(end)])
            xticks(unique(round(linspace(xVals(1),xVals(end),max(ceil(length(xVals)/8),5)),2)))
            
            xlabel('Z-Coordinate (mm)','Interpreter','latex')
            ylabel('Time (ps)','Interpreter','latex')
            
        elseif pInd{1} == 2
            
            xVals = (pInd{2}/max(pInd{2}))*this.props_len/10^-3;
            tVals = this.list_lambda(pInd{3},pInd{4})/10^-9;
            
            surf( xVals, tVals,...
                abs(squeeze( this.eField(pInd{:}) )').^2 /...
                max(max(abs(squeeze( this.eField(pInd{:}) )').^2)),...
                'EdgeColor','none','FaceLighting','gouraud');
%             view([45 30])
             view([0 90])
            
            xlim([xVals(1) xVals(end)])
            ylim([tVals(end) tVals(1)])
            xticks(unique(round(linspace(xVals(1),xVals(end),max(ceil(length(xVals)/8),5)),2)))
            
            xlabel('Z-Coordinate (mm)','Interpreter','latex')
            ylabel('Wavelength (nm)','Interpreter','latex')
                        
        end
        
        zlabel('Intesity (Norm.)','Interpreter','latex')
        
        
    end
    
end

set(gca,'FontSize',20);

end


function outs = parsePlotType(fieldInd,fieldNums,Nz)

outs = cell(1,3);

if ischar(fieldInd{1})
    if fieldInd{1} == 't'
        outs{1} = 1;
    elseif fieldInd{1} == 'w'
        outs{1} = 2;
    elseif strcmpi(fieldInd{1},'1') || strcmpi(fieldInd{1},'2')
        outs{1} = fieldInd{1};
    else
        error('Wrong first field type specified. Either ''t'' or ''w''.')
    end
else
    outs{1} = fieldInd{1};
    if ~(outs{1} == 1 || outs{1} == 2)
        error('Wrong first field type specified. Either 1 or 2.')
    end
end

if isa(fieldInd{2},'double')
    outs{2} = fieldInd{2};
elseif isa(fieldInd{2},'char')
    outs{2} = str2num(fieldInd{2}); %#ok<*ST2NM>
end

if any(outs{2} <= 0) || any(outs{2} > Nz+1)
    error(['Z-Step out of range. Must be between 1 and ',num2str(Nz+1)])
end


outs{3} = fieldNums;


end

function outs = findRange(this,fieldType,fieldRange)

if strcmpi(fieldRange,'all')
    
    fieldRange = [-inf,inf];
    
else
    
    if isa(fieldRange,'cell')
        fieldRange = parseUnits(fieldRange);
        fieldRange = cell2mat(fieldRange);
    end
    
    if fieldType{1} == 1
        if floor(log10(abs(fieldRange(1)))) >= 0
            fieldRange = fieldRange*10^-12;
        end
    elseif fieldType{1} == 2
        if floor(log10(abs(fieldRange(1)))) >= 0
            fieldRange = fieldRange*10^-9;
        end
    end
    
end

if fieldType{1} == 1
    
    outs{1} = find(...
        this.list_t >= fieldRange(1) &...
        this.list_t <= fieldRange(2) ...
        );
    
else
    
    outs{1} = find(...
        this.list_lambda(fieldType{3},:) >= fieldRange(1) &...
        this.list_lambda(fieldType{3},:) <= fieldRange(2) ...
        );
    
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

end

function checkInputs(this,fieldInd,fieldNum,fieldRange,norm)

if ~isa(fieldInd,'cell') || length(fieldInd) ~= 2
    error('Field type must be a 2x1 cell')
end

if length(fieldNum) > 1
    error('To many fields specified');
end
if ~(fieldNum <= size(this.eField,3) && fieldNum >= 1)
    error('Field number must evaluate to real fields');
end

if length(fieldRange) ~= 2 && ~strcmpi(fieldRange,'all')
    error('Plot range must contain exactly 2 points')
end

if ~(norm == 0 || norm == 1)
    error('Norm input is a boolean valued input: 0 or 1.')
end

end
