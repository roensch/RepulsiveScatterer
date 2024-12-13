function F = set_parameters(F,p,default)

fnames =  fieldnames(F);
for counter = 1:length(fnames);
    fname = fnames{counter};
    if isfield(p,fname)
        x=getfield(p,fname);
        if isstruct(x) & ~isempty(getfield(F,fname))
            x = complete_parameter_structure(x,getfield(F,fname));
        end
        F = setfield(F,fname,x);
    end
end
fnames = fieldnames(p);
for counter = 1:length(fnames)
    fname =  fnames{counter};
    if ~any(strcmp(properties(F),fname))
        fprintf(['WARNING: ' fname ' is defined in the parameter structure, but it is not a property of the class!\n']);
    end
end
if nargin > 2
    fnames =  fieldnames(default);
    for counter = 1:length(fnames);
        fname = fnames{counter};
        if isfield(p,fname)
            x=getfield(p,fname);
            if isstruct(x)
              x = complete_parameter_structure(x,getfield(default,fname));
            end
            F = setfield(F,fname,x);
        else
            x=getfield(default,fname);
            F = setfield(F,fname,x);
        end
    end
end
end