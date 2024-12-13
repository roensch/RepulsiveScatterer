classdef regmethod
    %% base class for regularization methods
    properties
        method = 'regmethod';
        % the higher the value, the more messages are given
        verbose = 1;
        % iteration step
        it_step = 0;
        % numbers of the reconstructions to be plotted
        plot_steps = [-1 0:1:100]
        % numbers of reconstructions to be stored
        store_rec = [1:100]
        stoprule_name = 'discrepancy';
        stoprule_par =struct('N_max_it',50,'tau',2);
        stoprule;
    end
    methods
        function R = regmethod(par)
            R = set_parameters(R,par);
            R.stoprule = feval(R.stoprule_name,R.stoprule_par);
        end
        
        function [stat,F] = error_stat_output(R,F,stat,x_k,y_obs,y_k,x_0)
            if R.it_step==0
                stat.Xerr=[];
                stat.Yerr=[];
                if ~isempty(R.store_rec)
                    stat.x = {};
                end
            end
            stat.Yerr = [stat.Yerr,sqrt(F.dot_Y(y_obs-y_k, y_obs-y_k))];            
            R.condPrint(1,'It.%i: resi=%1.3e\n',R.it_step,stat.Yerr(end));
            R.condPrint(1,'----------------------------------------------\n');
            if ~isempty(intersect(R.store_rec,R.it_step))
                stat.x = [stat.x,x_k];
            end
        end
        
        function condPrint(R,verbose_nr,varargin)
            if R.verbose>=verbose_nr
                if isnumeric(varargin{1})
                    fprintf(varargin{1},varargin{2:end});
                    if varargin{1}>=3
                        fprintf(varargin{2:end});
                    end
                else
                    fprintf(varargin{:})
                end
            end
        end
    end
end
