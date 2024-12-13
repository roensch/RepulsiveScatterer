classdef discrepancy
    % Implements Morozov's discrepancy principle
    
    properties
        ind = -1;
        rec = [];
        tau = 2;
        N_max_it = inf;
        name = 'discrepancy principle';
    end
    
    methods
        function stoprule = discrepancy(par)
            stoprule = set_parameters(stoprule,par);
        end
        
        function [muststop,stoprule] = stop(stoprule, F,xn,yn,data,R)
            stoprule.ind = stoprule.ind+1;
            stoprule.rec = xn;
            discr = sqrt(F.dot_Y(data-yn, data-yn));
            
            if (R.it_step < stoprule.N_max_it)
                if (discr > stoprule.tau*F.noiselevel)
                    muststop = false;
                    R.condPrint(2,'Discrepancy criterion failed as %f > %f.\n',discr, stoprule.tau*F.noiselevel);
                else
                    muststop = true;
                    R.condPrint(1,'Iteration stopped as discrepancy criterion is fulfilled (%f <= %f).\n',discr, stoprule.tau*F.noiselevel);
                end;
            else
                muststop = true;
                R.condPrint(1,'Iteration stopped as maximum number of iterations reached.\n')
            end
        end
        
        function [ind,res] = select_index(stoprule,F)
            ind = stoprule.ind;
            res = stoprule.rec;
        end
    end
end