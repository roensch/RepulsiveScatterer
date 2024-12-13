function p_new = complete_parameter_structure(p,p_ref)
% compares parameters specified in a structure p to a reference structure
% p_ref and returns a parameter structure p_new in which missing entries in p
% are taken from p_ref

p_new = p_ref;
fnames =  fieldnames(p_ref);
for counter = 1:length(fnames);
    fname = fnames{counter};
    if isfield(p,fname)
        x=getfield(p,fname);
        if ~isstruct(x)
            p_new = setfield(p_new,fname,x);        
        else  % recursive application for structs
            p_new = setfield(p_new,fname, ...
                 complete_parameter_structure(x,getfield(p_ref,fname)));
        end
    end
end
fnames = fieldnames(p);
for counter = 1:length(fnames)
    fname =  fnames{counter};
    if ~isfield(p_ref,fname)
        fprintf(['WARNING: No reference value specified for parameter ' fname '!\n']);
    end
end
end