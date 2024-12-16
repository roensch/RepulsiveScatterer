
function [prbl,reg] = reconstruction_parameters(experiment,save_flag)

    %%  parameters of the operator
    prbl.op_name = 'DirichletTri';
    
    prbl.metric = 'Tangent_Point';
    
    reg.max_remesher_discrepancy = 1;
    
    reg.Armijo_constant = 1/2;
    
    switch experiment
        case 'affinely_deformed_sphere'
            % create the true surface
            A = [0.19,0.3917,0.05;0.574,0.321,0.313;0.262,0.27,0.661]*[1.5,0,0;0,1.5,0;0,0,1/2];
            data_true = load("./Meshes/Sphere_00040560.mat").cell;
            Vcoord = data_true{1}*A;
            Conn = double(data_true{2});
            if det(A) < 0
                Conn = Conn(:,[1,3,2]);
            else
            end
            T_true = triangulation(Conn,Vcoord);
            data = load("./Meshes/Sphere_00040560.mat").cell;
            Vcoord = (1/2)*Vcoord;
            Conn = double(data{2});
            T = triangulation(Conn,Vcoord);
            prbl.init_guess = T;
    
            % set the iteration parameters
            prbl.noiselevel = 0.01;
            
            reg.alpha_0 = 1e-3;
            
            prbl.kappa = [pi,2*pi]; 
    
        case 'affinely_deformed_torus'
            % create the true surface
            A =[0.8420,    0.5763 ,  0.3813;...
                0.2817  ,  0.0718   , 0.6822;...
                0.1482   , 0.9724  ,  0.9148]*[1,0,0; 0,0,-1;0,1,0]*diag([1,1,1/2]);
            data_true = load("./Meshes/Torus_00038400.mat").Expression1;
            C = data_true{1};
    
            Vcoord = C * A;
            Conn = double(data_true{2});
            if det(A) < 0
                Conn = Conn(:,[1,3,2]);
            end
            T_true = triangulation(Conn,Vcoord);
    
            prbl.true_surface = T_true;
    
            % create the initial guess
            data = load("./Meshes/Torus_00038400.mat").Expression1;
            Vcoord = (1/2) * data{1};
            Conn = double(data{2});
            T = triangulation(Conn,Vcoord);
            prbl.init_guess = T;
            
            % set the iteration parameters
            prbl.noiselevel = 0.01;
            
            reg.alpha_0 = 1e-3;
            
            prbl.kappa = [pi,2*pi];  
    
        case 'snake_torus'
            % create the true surface
            T_true = stlread('./Meshes/snake_torus.stl');
            T_true = triangulation(T_true.ConnectivityList,(1/2)*T_true.Points);
            prbl.true_surface = T_true;
    
            % create the initial guess
            data = load("./Meshes/Torus_00038400.mat").Expression1;
            Vcoord = (1/8)*data{1};
            Conn = double(data{2});
            T = triangulation(Conn,Vcoord);
            prbl.init_guess = T;
    
            % set the iteration parameters
            prbl.noiselevel = 0.01;
            
            reg.alpha_0 = 1e-3;
            
            prbl.kappa = [pi,2*pi];  
    
        case 'bob'
            % create the true surface
            turn = [1,0,0;...
                0,0,1;...
                0,1,0];
            data_true = load("./Meshes/Bob_00171008.mat").Expression1;
            Vcoord = data_true{1} * turn;
            Conn = double(data_true{2})* turn;
    
            T_true = triangulation(Conn,Vcoord);
            prbl.true_surface = T_true;
    
            % create the initial guess
            data = load("./Meshes/Torus_00038400.mat").Expression1;
            Vcoord = data{1}*[1/2,0,0; 0,1/2,0;0,0,1];
            Conn = double(data{2});
            T = triangulation(Conn,Vcoord);
            prbl.init_guess = T;
            
            % set the iteration parameters
            prbl.noiselevel = 0.01;
            
            reg.alpha_0 = 1e-3;
            
            prbl.kappa = [pi,2*pi]; 
    
        case 'bob_sphere'
            % create the true surface
            turn = [1,0,0;...
                0,0,1;...
                0,1,0];
            data_true = load("./Meshes/Bob_00171008.mat").Expression1;
            Vcoord = (1/2) * data_true{1} * turn;
            Conn = double(data_true{2})* turn;
    
            T_true = triangulation(Conn,Vcoord);
            prbl.true_surface = T_true;
    
            % create the initial guess
            data = load("./Meshes/Sphere_00081920.mat").Expression1;
            Vcoord = data{1};
            Conn = double(data{2});
            T = triangulation(Conn,Vcoord);
            prbl.init_guess = T;
            
            % set the iteration parameters
            prbl.noiselevel = 0.01;
            
            reg.alpha_0 = 1e-3;
            
            prbl.kappa = [pi,2*pi];
    
        case 'spot'
            turn = [1,0,0;...
                0,0,1;...
                0,1,0];
    
            % create the true surface
            data_true = load("./Meshes/Spot_00093696.mat").Expression1;
            Vcoord = data_true{1} * turn;
            Conn = double(data_true{2})* turn;
            T_true = triangulation(Conn,Vcoord);
            prbl.true_surface = T_true;
    
            % create the initial guess
            data_true = load("./Meshes/Sphere_00081920.mat").Expression1;
            Vcoord = data_true{1};
            Conn = double(data_true{2});
            T = triangulation(Conn,Vcoord);
            prbl.init_guess = T;
    
            % set the iteration parameters
            prbl.noiselevel = 0.01;
            
            reg.alpha_0 = 1e-3;
            
            prbl.kappa = [pi,2*pi];
    
         case 'blub'
            turn = [1,0,0;...
                0,0,1;...
                0,1,0];
    
            % create the true surface
            data_true = load("./Meshes/Blub_00227328.mat").Expression1;
            Vcoord = data_true{1} * turn;
            Conn = double(data_true{2})* turn;
            T_true = triangulation(Conn,Vcoord);
            prbl.true_surface = T_true;
    
            % create the initial guess
            data_true = load("./Meshes/Sphere_00081920.mat").Expression1;
            Vcoord = data_true{1};
            Conn = double(data_true{2});
            T = triangulation(Conn,Vcoord);
            prbl.init_guess = T;
    
            % set the iteration parameters
            prbl.noiselevel = 0.01;
            
            reg.alpha_0 = 1e-3;
            
            prbl.kappa = [pi,2*pi];
    
        case 'bunny'
            % create the true surface
            T_t = stlread('./Meshes/bunny.stl');
            Vcoord = (1/50) * T_t.Points - repmat([0,0,1/2],size(T_t.Points,1),1);
            T_true = triangulation(T_t.ConnectivityList,Vcoord);
            prbl.true_surface = T_true;
    
    %         create the initial guess
            data_true = load("./Meshes/Sphere_00040560.mat").cell;
            Vcoord = data_true{1};
            Conn = double(data_true{2});
    
            T = triangulation(Conn,Vcoord);
            prbl.init_guess = T;
    
            % set the iteration parameters
            prbl.noiselevel = 0.01;
            
            reg.alpha_0 = 1e-3;
            
            prbl.kappa = [pi,2*pi];
            
        case 'triceratops'
            % create the true surface
            data_true = load("./Meshes/Triceratops_00090560.mat").Expression1;
            Vcoord = (1/4)*data_true{1};
            Conn = double(data_true{2});
            T_true = triangulation(Conn,Vcoord);
            prbl.true_surface = T_true;
    
            % create the initial guess
            data_true = load("./Meshes/Sphere_00081920.mat").Expression1;
            Vcoord = data_true{1};
            Conn = double(data_true{2});
            T = triangulation(Conn,Vcoord);
            prbl.init_guess = T;
    
            % set the iteration parameters
            prbl.noiselevel = 0.01;
            
            reg.alpha_0 = 1e-3;
            
            prbl.kappa = [pi,2*pi];
    
        case 'dolphin'
            % create the true surface
            T_t = stlread('./Meshes/dolphin_remeshed_2.stl');
            Vcoord = T_t.Points;
            T_true = triangulation(T_t.ConnectivityList,Vcoord);
            prbl.true_surface = T_true;
    %         trisurf(T_t.ConnectivityList,Vcoord(:,1),Vcoord(:,2),Vcoord(:,3))
    
    %         create the initial guess
            data_true = load("./Meshes/Sphere_00081920.mat").Expression1;
            Vcoord = data_true{1};
            Conn = double(data_true{2});
    
            T = triangulation(Conn,Vcoord);
            prbl.init_guess = T;
    
            % set the iteration parameters
            prbl.noiselevel = 0.01;
            
            reg.alpha_0 = 1e-3;
            
            prbl.kappa = [pi,2*pi];
            
    end
    
    I = readmatrix("./Meshes/FeketePoints_00000016V_reduced.txt");
    prbl.inc_dir = I;
    
    %%  parameters of the regularization method
    % regularization method
    reg.method = 'GaussNewton';
    
    % set maximal iterations and discrepancy parameter
    reg.stoprule_par = struct('N_max_it',100,'tau',3);
    reg.verbose = 2;
    
    %% Set saving path if requested
    reg.save_flag = save_flag;

    if ( (save_flag ~= "") && not(isfolder('./Reconstructions')) )
        mkdir('./Reconstructions');
    end
end
