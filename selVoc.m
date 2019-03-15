function DataSel = selVoc(Res)
%% Get the data ready
% Need to get rid of mlnoise sections and whine sections when they
% exist. I construct a vector of indices of the right sections
DataSel=zeros(1,length(Res.VocType));
nvoc=0;
voctype=Res.VocType;
for ii=1:length(Res.VocType)
    if strcmp(voctype{ii}, 'Ag')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'Be')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'DC')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'Di')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'LT')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'Ne')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'Te')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'Th')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
    elseif strcmp(voctype{ii}, 'song')
        nvoc=nvoc+1;
        DataSel(nvoc)=ii;
        %     elseif strcmp(voctype{dd}, 'Wh')
        %         nvoc=nvoc+1;
        %         DataSel(nvoc)=dd;
    end
end
DataSel=DataSel(1:nvoc);
return 