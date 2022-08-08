function Results = aeolianBedformSurfaceModel()
%A surface model for aeolian dune topography
% Mathematical Geosciences
% Travis Swanson 1*, David Mohrig 1, Gary Kocurek 1, Man Liang 1
% 1 Department of Geosciences, Jackson School of Geosciences, 
% The University of Texas at Austin, 2305 Speedway, Stop C1160, 
% Austin, TX 78712-1692
% *travis.swanson@utexas.edu
%
%
% Results structure has fields:
%   etas            3D array of topography, 3rd dimension is time
%   name            simulation name (string)
%   meanHeight      Time series of mean bedform height 
%   stdHeight       Time series of bedform height standard deviation
%   meanWavelength  Time series of mean wavelength
%   stdWavelength   Time series of bedform wavelength standard deviation
%   crestLogical    Logical array of crest locations
%   iAngle          3D array surface dip direction
%   dipAngle        3D array surface dip angle
%   terminations    Time series of number of crest terminations
%
% run >>help post 
% for more information


%clear out memory
clear variables

fprintf('Cases from the manuscript are in order below.\n');
fprintf('1) uniper - unimodal transport regime, periodic boundary condtions\n');
fprintf('2) biper - bimodal transport regime, periodic boundary condtions\n');
fprintf('3) uniline - unimodal transport regime, line source area\n');
fprintf('4) biline - bimodal transport regime, line source area\n');
inputSimNum = input('please input simulation number (1-4) : ');

%choose a simulation to run (1 through 4).
simulation_number = inputSimNum;

%create input for simulation
fprintf('model input \n');
Input = input_file(simulation_number);
%create 3D array of aeolian bedform topography
fprintf('running aeolian surface model \n');
etas = ae_interface(Input);
%post processing
fprintf('post processing model results \n');
Results = post(etas,Input);

end

%% Model input creation

function Input = input_file(sim_i)


    %Each set of input conditions are pre-created for each simulation (1 through 4).
    switch sim_i
        case  1
                     Input.simulation = 'periodic unimodal';
        case  2
                     Input.simulation = 'periodic bimodal';
        case  3
                     Input.simulation = 'line unimodal';
        case  4
                     Input.simulation = 'line bimodal';
    end


    switch Input.simulation

        case 'line unimodal';
            Input.boundaryCondition = 'Flux';
            sedSource = 'side';
            sedInitialCondition = 'Flat';
            windMode = 'Unimodal';
            Input.startRecording = 10000; %time step to start recording
            Input.nt=100000; % number of time steps to run
            Input.DetermOrStoch = 0; % 0 = constant shear stress threshold 1 = stochastic shear stress threshold

        case 'line bimodal';
            Input.boundaryCondition = 'Flux';
            sedSource = 'side';
            sedInitialCondition = 'Flat';
            windMode = 'Bimodal';
            Input.startRecording = 10000;
            Input.nt=100000; % number of time steps to run
            Input.DetermOrStoch = 0;

        case 'periodic unimodal';
            Input.boundaryCondition = 'Periodic';
            sedSource = 'side';
            sedInitialCondition = 'Area';
            windMode = 'Unimodal';
            Input.startRecording = 0;
            Input.nt=20000; % number of time steps to run
            Input.DetermOrStoch = 0;

        case 'periodic bimodal';

            Input.boundaryCondition = 'Periodic';
            sedSource = 'side';
            sedInitialCondition = 'Area';
            windMode = 'Bimodal';
            Input.startRecording = 0;
            Input.nt= 30000; % number of time steps to run
            Input.DetermOrStoch = 0;
    end


    Input.record_eta = true; %record surfaces?
    Input.SampleFreq = 10; %sample every __ time steps (saves memory)



    %turn on or off first order term of Exner
    Input.oot = 1;
    %make the grid
    l_x = 1000; %length of Easting axis
    Input.dx = 10; %spatial step (dx or dy)
    x = 0:Input.dx:l_x; % length of domain (length of 1D profiles)
    y = 0:Input.dx:2*l_x; %Northing axis is twice the easting axis
    Input.x = x;
    Input.y = y;

    Input.dt = 0.7; %time step, dt

    %create 2D arrays for plotting surface.
    [Input.X,Input.Y] = meshgrid(y,x);



    %transport conditions
    %init program varibles

    Input.A =0.1; %shape parameter 1
    Input.B = 3; %shape parameter 2
    Input.E = 10; %avalanche fraction (1= avalanche occurs in one timestep)

    if Input.DetermOrStoch == 1
        %this sets the shear stress threshold @ 0
        Input.tcrit = 0;
    else
        %this makes a stochastic shear stress threshold (each node has a
        %randomly chosen value from this function).
        Input.tcrit = @(x) 3*rand(size(x));
    end

    Input.m = 1; % sediment transport (power law) coeff.
    Input.n = 1.5; % sediment transport (power law) exp.
    % Input.eta = 1*rand(size(Input.X)); %initial condition
    Input.tc = tand(32); % angle of repose
    Input.p = 0.4; % porosity
    Input.Dg = 0.08;

    %initial conditions
    switch sedInitialCondition
        case 'Area'
            Input.eta = 1*rand(size(Input.X));
        case 'Flat'
            Input.eta = zeros(size(Input.X));
    end

    if Input.record_eta == 1
        Input.etas=zeros(size(Input.eta,1),size(Input.eta,2),(Input.nt - Input.startRecording)/Input.SampleFreq);
        Input.cetas=zeros(size(Input.eta,1),size(Input.eta,2),(Input.nt - Input.startRecording)/Input.SampleFreq);
    end

    switch sedSource
        case 'side'
          fbc = @(x) 10*rand(size(Input.eta,1),2);
          Input.fbc = fbc;
    end


    %Create a transport direction for each time step of simulation.

    switch windMode
        case 'Unimodal'
               %This is the unimodal transport regime
               % mean = 0, std = 15
               windDirection = 0 + 15*randn(1,Input.nt);
        case 'Bimodal'
    %         windDirection = zeros(1,Input.nt);
            windDirection = [];
            w_interval = 20; %number of time steps for each annual wind
            wd = 1;
            for i = 1:w_interval:Input.nt
                windDirection = [windDirection wd*80 + 20*randn(1,w_interval)];
                wd = -wd;
            end
             %bimodal_wind = [60 + 20*randn(1,ceil(Input.nt/2)),  - 60 + 20*randn(1,ceil(Input.nt/2))];
             %windDirection = bimodal_wind(randperm(length(bimodal_wind)));
        case 'Unidirectional'
            windDirection = 0*ones(1,Input.nt);
        case 'Bidirectional'
            windDirection = repmat([60,-60],[1 Input.nt/2]); 
    end

    %create wind regime
    ws = 1; % this COULD be made time variable as well; simply make this vector of length Input.nt
    Input.tbx = ws*cosd(windDirection);
    Input.tby = ws*sind(windDirection);


    % create indexing vectors for finite differences
    switch Input.boundaryCondition

        case 'Periodic'
            Input.x_index = 1:size(Input.eta,2);
            Input.y_index = 1:size(Input.eta,1);
            Input.x_p1_index = [2:size(Input.eta,2) 1];
            Input.x_m1_index = [size(Input.eta,2) 1:size(Input.eta,2)-1];
            Input.y_p1_index = [2:size(Input.eta,1) 1];
            Input.y_m1_index = [size(Input.eta,1) 1:(size(Input.eta,1)-1)];
            Input.fluxBoundary = false;

         case 'Walls'
            Input.x_index = 1:length(Input.eta);
            Input.y_index = 1:length(Input.eta);
            Input.x_p1_index = [2:length(Input.eta) 1];
            Input.x_m1_index = [length(Input.eta) 1:length(Input.eta)-1];
            Input.y_p1_index = [2:length(Input.eta) length(Input.eta)];
            Input.y_m1_index = [1 1:(length(Input.eta)-1)];
            Input.fluxBoundary = false;

         case 'Flux'
             % The outter ring sets the transport conditions inward.
             %full eta n x n
            Input.x_index = 2:length(Input.eta)-1;
            Input.y_index = 2:size(Input.eta,1)-1;
            Input.x_p1_index = 3:length(Input.eta);
            Input.x_m1_index = 1:length(Input.eta)-2;
            Input.y_p1_index = 3:size(Input.eta,1);
            Input.y_m1_index = 1:(size(Input.eta,1)-2); 
            %n-1 x n-1
            %n-2 x n-2
            Input.fluxBoundary = true;

    end
end

%% Aeolian bedform surface model

  function [etas] = ae_interface(Inputs)
% Model for generation of an aeolian dune field surface from noise, with
% periodic/flux boundary conditions. 

%place the input structure fields into local variables

createLocalVar(Inputs);

ic = 1;


if fluxBoundary == 0
    for j = 1:nt

         if tbx(j) > 0 %This is the condition for upwind in the X direction
                  tx = tbx(j).*(oot + A.*(eta(y_index,x_index))+B.*((eta(y_index,x_index)-eta(y_index,x_m1_index))./dx));
             elseif tbx(j) < 0%This is the condition for down wind in the X direction
                    tx = -tbx(j).*(oot +A.*(eta(y_index,x_index))+B.*((eta(y_index,x_index)-eta(y_index,x_p1_index))./dx));
             elseif tbx(j) == 0
                 tx = zeros(size(eta(y_index,x_index)));
         end
    % 
        if tby(j) > 0 %This is the condition for upwind in the X direction
                  ty = tby(j).*(oot +A.*(eta(y_index,x_index))+B.*((eta(y_index,x_index)-eta(y_m1_index,x_index))./dx));
            elseif tby(j) < 0 %This is the condition for down wind in the X direction
                  ty = -tby(j).*(oot +A.*(eta(y_index,x_index))+B.*((eta(y_index,x_index)-eta(y_p1_index,x_index))./dx));
            elseif tby(j) == 0
                ty = zeros(size(eta(y_index,x_index)));
        end

        if DetermOrStoch == 1 %constant threshold
             tx(tx<tcrit) = 0;
             ty(ty<tcrit) = 0;
        else %stochastic threshold
             tx(tx<tcrit(tx)) = 0;
             ty(ty<tcrit(ty)) = 0;
        end

       %let's compute a crieteria for -+X and +- Y directions.
       %and relax _any_ slope that becomes greater than threshold.

        if tbx(j) < 0 % both down wind
                  avhCrit_dwx = (eta(y_index,x_index)-eta(y_index,x_m1_index))./dx;
                  qa_x =  E*(avhCrit_dwx.^2 - tc^2).*avhCrit_dwx.*((avhCrit_dwx) > tc);
            elseif tbx(j) > 0  %upwind x downwind y
                avhCrit_uwx = (eta(y_index,x_index)-eta(y_index,x_p1_index))./dx;
                qa_x = E*(avhCrit_uwx.^2 - tc^2).*avhCrit_uwx.*((avhCrit_uwx) > tc);
            elseif tbx(j) == 0
                qa_x = zeros(size(tx));
        end

         if tby(j) > 0 % downwind x upwind y
                 avhCrit_uwy = (eta(y_index,x_index)-eta(y_p1_index,x_index))./dx;
                 qa_y = E*(avhCrit_uwy.^2 - tc^2).*avhCrit_uwy.*(avhCrit_uwy > tc);
             elseif tby(j) < 0 %both upwind
                 avhCrit_dwy = (eta(y_index,x_index)-eta(y_m1_index,x_index))./dx;
                 qa_y = E*(avhCrit_dwy.^2 - tc^2).*avhCrit_dwy.*(avhCrit_dwy > tc);
              elseif tby(j) == 0
                  qa_y = zeros(size(ty));
         end

            qxc = sqrt((m.*tx.^n + qa_x ).^2 + (m.*ty.^n  + qa_y).^2);

         if tbx(j) < 0 % both down wind
                 xcomp = qxc(y_index,x_index)-qxc(y_index,x_p1_index);
            elseif tbx(j) > 0  %upwind x downwind y
                xcomp = qxc(y_index,x_index)-qxc(y_index,x_m1_index);
             elseif tbx(j) == 0
                 xcomp = 0;
         end

         if tby(j) > 0 % downwind x upwind y
                ycomp = qxc(y_index,x_index)-qxc(y_m1_index,x_index);
             elseif tby(j) < 0 %both upwind
                ycomp =  qxc(y_index,x_index)-qxc(y_p1_index,x_index);
             elseif tby(j) == 0
                ycomp = 0;
         end

        a_eta = (-dt/((1-p)*dx)).*abs(tbx(j)).*xcomp  + (-dt/((1-p)*dx)).*abs(tby(j)).*ycomp;

        %topographic diffusion 
        d_eta = (dt*Dg/(dx^2)).*(eta(y_index,x_p1_index)+eta(y_index,x_m1_index) + eta(y_m1_index,x_index) + eta(y_p1_index,x_index)  - 4.*eta(y_index,x_index));
        ceta = d_eta + a_eta;
        eta(y_index,x_index)  = eta(y_index,x_index) +ceta;

        if record_eta == 1 && mod(j,SampleFreq) == 0
                etas(:,:,ic)=eta;
                ic = ic + 1;
        end

        if mod(j,10000) == 0 %display surface every 10^4 time steps.
            surf(X(y_index,x_index),Y(y_index,x_index),eta(y_index,x_index)); axis equal;  view(2); shading interp; colormap jet; drawnow;
            title('Surface Elevation');
            drawnow
        end
    end

else %case where there is a flux boundary condition %%%%%%%%%%%%%%%%%%%%%%
    
    for j = 1:nt
    

            eta(:,1:2) = fbc(eta);
            eta(1:2,:) = repmat(eta(3,:),2,1);
            eta(end-1:end,:) = repmat(eta(end-2,:),2,1);

        
       %shear stress calculations need to access all nodes, nxn, gives
       %results for n-1 x n-1

         if tbx(j) > 0 %This is the condition for upwind in the X direction
                  tx = tbx(j).*(oot + A.*(eta(y_index,x_index))+B.*((eta(y_index,x_index)-eta(y_index,x_m1_index))./dx));
             elseif tbx(j) < 0%This is the condition for down wind in the X direction
                    tx = -tbx(j).*(oot +A.*(eta(y_index,x_index))+B.*((eta(y_index,x_index)-eta(y_index,x_p1_index))./dx));
             elseif tbx(j) == 0
                 tx = zeros(size(eta(y_index,x_index)));
         end
    % 
        if tby(j) > 0 %This is the condition for upwind in the X direction
                  ty = tby(j).*(oot +A.*(eta(y_index,x_index))+B.*((eta(y_index,x_index)-eta(y_m1_index,x_index))./dx));
            elseif tby(j) < 0 %This is the condition for down wind in the X direction
                  ty = -tby(j).*(oot +A.*(eta(y_index,x_index))+B.*((eta(y_index,x_index)-eta(y_p1_index,x_index))./dx));
            elseif tby(j) == 0
                ty = zeros(size(eta(y_index,x_index)));
        end


        if DetermOrStoch == 1 %constant threshold
             tx(tx<tcrit) = 0;
             ty(ty<tcrit) = 0;
        else %stochastic threshold
             tx(tx<tcrit(tx)) = 0;
             ty(ty<tcrit(ty)) = 0;
        end


       %let's compute a crieteria for -+X and +- Y directions.
       %and relax _any_ slope that becomes greater than threshold.
      

        if tbx(j) < 0 % both down wind
              avhCrit_dwx = (eta(y_index,x_index)-eta(y_index,x_m1_index))./dx;
              qa_x =  (avhCrit_dwx.^2 - tc^2).*avhCrit_dwx.*((avhCrit_dwx) > tc);
        elseif tbx(j) > 0  %upwind x downwind y
            avhCrit_uwx = (eta(y_index,x_index)-eta(y_index,x_p1_index))./dx;
            qa_x = (avhCrit_uwx.^2 - tc^2).*avhCrit_uwx.*((avhCrit_uwx) > tc);
        elseif tbx(j) == 0
            qa_x = zeros(size(tx));
        end
    %     
         if tby(j) > 0 % downwind x upwind y
             avhCrit_uwy = (eta(y_index,x_index)-eta(y_p1_index,x_index))./dx;
             qa_y = (avhCrit_uwy.^2 - tc^2).*avhCrit_uwy.*(avhCrit_uwy > tc);
         elseif tby(j) < 0 %both upwind
             avhCrit_dwy = (eta(y_index,x_index)-eta(y_m1_index,x_index))./dx;
             qa_y = (avhCrit_dwy.^2 - tc^2).*avhCrit_dwy.*(avhCrit_dwy > tc);
          elseif tby(j) == 0
              qa_y = zeros(size(ty));
         end
         
         %qxc is n-1 x n-1

        qxc = sqrt((m.*tx.^n + qa_x ).^2 + (m.*ty.^n  + qa_y).^2);
  
        %xcomp, ycomp will be n-2 x n-2; but need to access n-1 x n-1; the
        %dimensions of qxc

         if tbx(j) < 0 % both down wind
                 xcomp = qxc(:,1:end-1)-qxc(:,2:end);
            elseif tbx(j) > 0  %upwind x downwind y
%                 xcomp = qxc(y_index,x_index(2:end))-qxc(y_index,x_m1_index(2:end));
                xcomp = qxc(:,2:end)-qxc(:,1:end-1);
             elseif tbx(j) == 0
                 xcomp = 0;
         end

         if tby(j) > 0 % downwind x upwind y
                ycomp = qxc(2:end,:)-qxc(1:end-1,:);
             elseif tby(j) < 0 %both upwind
                ycomp =  qxc(1:end-1,:)-qxc(2:end,:);
             elseif tby(j) == 0
                ycomp = 0;
         end
         
        %topographic diffusion
        d_eta = (dt*Dg/(dx^2)).*(eta(y_index,x_p1_index)+eta(y_index,x_m1_index) + eta(y_m1_index,x_index) + eta(y_p1_index,x_index)  - 4.*eta(y_index,x_index));
       
       %depending on the transport vector; choose your upwind/downwind yes.
       if tby(j) == 0 && tbx(j) == 0
               a_eta = 0;
           elseif tby(j) > 0 && tbx(j) == 0 % y+
                a_eta = (-dt/((1-p)*dx)).*abs(tby(j)).*ycomp(1:end-1,2:end-1); 
           elseif tby(j) < 0 && tbx(j) == 0 % y-
                a_eta = (-dt/((1-p)*dx)).*abs(tby(j)).*ycomp(2:end,2:end-1); 
           elseif tby(j) == 0 && tbx(j) > 0 % x+
                a_eta = (-dt/((1-p)*dx)).*abs(tbx(j)).*xcomp(2:end-1,1:end-1); 
           elseif tby(j) == 0 && tbx(j) < 0 % x-
                a_eta = (-dt/((1-p)*dx)).*abs(tbx(j)).*xcomp(2:end-1,2:end);
           elseif tby(j) > 0 && tbx(j) > 0 % y+ x+
                a_eta = (-dt/((1-p)*dx)).*abs(tbx(j)).*xcomp(2:end-1,1:end-1)  + (-dt/((1-p)*dx)).*abs(tby(j)).*ycomp(1:end-1,2:end-1);
           elseif tby(j) > 0 && tbx(j) < 0 % y+ x-
                a_eta = (-dt/((1-p)*dx)).*abs(tbx(j)).*xcomp(2:end-1,2:end)  + (-dt/((1-p)*dx)).*abs(tby(j)).*ycomp(1:end-1,2:end-1);
           elseif tby(j) < 0 && tbx(j) < 0 % y- x-
                a_eta = (-dt/((1-p)*dx)).*abs(tbx(j)).*xcomp(2:end-1,2:end)  + (-dt/((1-p)*dx)).*abs(tby(j)).*ycomp(2:end,2:end-1);
           elseif tby(j) < 0 && tbx(j) > 0 % y- x+
                a_eta = (-dt/((1-p)*dx)).*abs(tbx(j)).*xcomp(2:end-1,1:end-1)  + (-dt/((1-p)*dx)).*abs(tby(j)).*ycomp(2:end,2:end-1); 
       end

         eta(3:end-2,3:end-2) = eta(3:end-2,3:end-2)  + d_eta(2:end-1,2:end-1) + a_eta;
        
        if record_eta == 1 && mod(j,SampleFreq) == 0 && j >= startRecording
            etas(:,:,ic)=eta;
            ic = ic + 1;
        end

        if mod(j,10000) == 0 %plot surface every 10^4 time steps
            surf(X(y_index,x_index),Y(y_index,x_index),eta(y_index,x_index)); axis equal;  view(2); shading interp; colormap jet; drawnow;     
            title('Surface Elevation');
            drawnow
        end
    end
end

if record_eta == 0 %return the last surface if told not to record.
    etas = eta;
end
  end
  
%% post-processing model results
 
function t = post(etas,Input)
 % Post processing function
 %
 % Input: etas  - 3D array of surface model output
 %        Input - structure created by input_file.m
 %
 % Output: 
 % t - Matlab structure with fields:
 %
 % t.etas           = (nxmxdt) surface model topography sampled from simulation
 %
 % t.name           = simulation name
 %
 % t.meanHeight     = mean values of bedform height for 1D transects of
 % simulation topography. This creates a 2d array (rows = time steps, 
 % colums = transect number)     
 %
 % t.stdHeight   = standard deviation values of bedform height for 1D transects of
 % simulation topography. This creates a 2d array (rows = time steps, 
 % colums = transect number)  
 %
 % t.meanTrough     = mean values of bedform trough elevation for 1D transects of
 % simulation topography. This creates a 2d array (rows = time steps, 
 % colums = transect number)   
 %
 % t.stdTrough   = standard deviation values of bedform trough elevation for 1D transects of
 % simulation topography. This creates a 2d array (rows = time steps, 
 % colums = transect number)  
 %
 % t.meanWavelength = mean values of bedform wavelength for 1D transects of
 % simulation topography. This creates a 2d array (rows = time steps, 
 % colums = transect number)  
 %
 % t.stdWavelength  = standard deviation values of bedform height for 1D transects of
 % simulation topography. This creates a 2d array (rows = time steps, 
 % colums = transect number) 
 %
 % t.crestLogical   = a 3D logical array. True entries are locations of a
 % dune crest in each time slice of t.etas
 %
 % t.iAngle         = 3D array of incidence angle (acute horizontal angle between
 % transport vector and bedform surface normal vector). each 2D slice is a
 % incidence angle map for each time slice in t.etas.
 %
 % t.dipAngle       = 3D array of surface dip magnitude. each 2D slice is a
 % dip angle map for each time slice in t.etas.
 %
 % t.terminations   = number of crest terminations for each time slice in t.etas
    
    if Input.fluxBoundary == false
        if strcmp(Input.simulation,'periodic bimodal')
            [mh,stdh,mt,stdt,mw,stdw,Li] = mean_props_through_time2d(reshape(etas,size(etas,2),size(etas,1),size(etas,3)));
        else
            [mh,stdh,mt,stdt,mw,stdw,Li] = mean_props_through_time2d(etas);
        end
        mh = mh-mt;
    
        [ina,da] = fast_surf2theta(Input.X,Input.Y,etas,[],Input.dx);

        endpts = TwoDCrestFind_logical_time(etas);
        
        t(1).etas = etas;
        t(1).name = Input.simulation;
        t(1).meanHeight = mh;
        t(1).stdHeight = stdh;
        t(1).meanTrough = mt;
        t(1).stdTrough = stdt;
        t(1).meanWavelength = mw;
        t(1).stdWavelength = stdw;
        t(1).crestLogical = Li;
        t(1).iAngle = ina;
        t(1).dipAngle = da;
        t(1).terminations = endpts;

    elseif Input.fluxBoundary == true;
    % plot flux boundary

        [mh,stdh,mw,stdw,Li] = mean_props_through_space2d(etas,Input.y,Input.dx,Input.simulation);

        [ina,da] = fast_surf2theta(Input.X,Input.Y,etas,[],Input.dx);
   
        endpts = TwoDCrestFind_logical_space(etas);
        endpts = mean(sum(endpts,3),1);
        
        t(1).etas = etas;
        t(1).name = Input.simulation;
        t(1).meanHeight = mh;
        t(1).stdHeight = stdh;
        t(1).meanWavelength = mw;
        t(1).stdWavelength = stdw;
        t(1).crestLogical = Li;
        t(1).iAngle = ina;
        t(1).dipAngle = da;
        t(1).terminations = endpts;

    end
 end
 
function [mh,stdh,mw,stdw,Li] = mean_props_through_space2d(etas,y,dx,simulation)
    % function to calculate mean values as a function of space for a 2D surface
    % model.
    %The input and output of this function are described in post.m
        Li = false(size(etas));
        H = zeros(size(etas));
        W = zeros(size(etas));

        if strcmp(simulation,'line bimodal')

            for i = 1:1:size(etas,3)
                for j = 1:1:size(etas,2)
                    [h,t,w,lipks,lipks_w] = OneDDunePropsSpace(squeeze(etas(:,j,i)));
                    H(:,j,i) = h-abs(t);
                    W(:,j,i) = w;

                    Li(:,j,i)=lipks; %logical array for binary image processing.  
                end
            end

            mh = (sum(sum(H,3),1))./sum(sum(H~=0,3),1);
            stdh = zeros(size(mh));
            for i = 1:length(mh)
                stdh(i) = sqrt(sum((reshape(squeeze(H(:,i,:)),1,[])-mh(i)).^2)./sum(reshape(H(:,i,:),1,[])~=0));
            end
            mw = (sum(sum(W,3),1))./sum(sum(W~=0,3),1);
            stdw = zeros(size(mh));
            for i = 1:length(mw)
                stdw(i) = sqrt(sum((reshape(squeeze(W(:,i,:)),1,[])-mw(i)).^2)./sum(reshape(W(:,i,:),1,[])~=0));
            end

        else %% the problem must be UNIMODAL

            for i = 1:1:size(etas,3)
                for j = 1:1:size(etas,1)
                    [h,t,w,lipks,lipks_w] = OneDDunePropsSpace(squeeze(etas(j,:,i)));
                    H(j,:,i) = h-abs(t);
                    W(j,:,i) = w;

                    Li(j,:,i)=lipks; %logical array for binary image processing.

                end
            end
            mh = (sum(sum(H,3),1))./sum(sum(H~=0,3),1);
            stdh = zeros(size(mh));
            for i = 1:length(mh)
                stdh(i) = sqrt(sum((reshape(squeeze(H(:,i,:)),1,[])-mh(i)).^2)./sum(reshape(H(:,i,:),1,[])~=0));
            end
            mw = (sum(sum(W,3),1))./sum(sum(W~=0,3),1);
            stdw = zeros(size(mh));
            for i = 1:length(mw)
                stdw(i) = sqrt(sum((reshape(squeeze(W(:,i,:)),1,[])-mw(i)).^2)./sum(reshape(W(:,i,:),1,[])~=0));
            end


        end
 end
 
function [mh,stdh,mt,stdt] = mean_props_through_time(etas)
    % function to calculate mean values as a function of time for a 2D surface
    % model. 1D only.
    %The input and output of this function are described in post.m

    mh = zeros(1,(size(etas,3)));
    mt = zeros(1,(size(etas,3)));
    stdh = zeros(1,(size(etas,3)));
    stdt = zeros(1,(size(etas,3)));
    ic = 1;

    for i = 1:1:size(etas,3)

            [mh(ic), stdh(ic), mt(ic), stdt(ic)] = OneDDuneProps(squeeze(etas(50,:,i)));
            ic = ic + 1;
    end
 end
 
function [mh,stdh,mt,stdt,mw,stdw,Li] = mean_props_through_time2d(etas)
    % function to calculate mean values as a function of time for a 2D surface
    % model.
    %The input and output of this function are described in post.m
    ic = 1;

    Li = false(size(etas));

    for i = 1:1:size(etas,3)
        for j = 1:1:size(etas,1)
            [mh(ic,j),stdh(ic,j),mt(ic,j),stdt(ic,j),mw(ic,j),stdw(ic,j),lipks] = OneDDuneProps(squeeze(etas(j,:,i)));

            Li(j,:,i)=lipks; %logical array for binary image processing.
        end
        ic = ic + 1;
    end
 end
 
function [mh,stdh,mt,stdt,mw,stdw,lipks] = OneDDuneProps(eta)
    %function to compute mean length scales of 1D dunes.
    %
    %MeanHeight is the mean amplitude of the train of dune crests
    %MeanWavelength is the peak-to-peak distance of the train
    %MeanStoss is the trough-to-peak distance in the direction of dune
    %migration
    %MeanLee is the trough-to-peak distance in the ~direction of dune
    %migration

    eta = eta + abs(min(eta));

    [~,Plocs] = findpeaks(eta);
    [~,Tlocs] = findpeaks(-eta);
    dx = 10;
    lipks = ismember(1:length(eta),Plocs);

    mh = mean(eta(Plocs));
    stdh = std(eta(Plocs));

    mt = mean(eta(Tlocs));
    stdt = std(eta(Tlocs));

    %calcualted by peaks
    mw = mean(diff(Plocs)*dx);
    stdw = std(diff(Plocs)*dx);

 end

function [htr,ttr,wtr,lipks,lipks_w] = OneDDunePropsSpace(eta)
    %function to compute mean length scales of 1D dunes.
    %
    %MeanHeight is the mean amplitude of the train of dune crests
    %MeanWavelength is the peak-to-peak distance of the train
    %MeanStoss is the trough-to-peak distance in the direction of dune
    %migration
    %MeanLee is the trough-to-peak distance in the ~direction of dune
    %migration
    htr = zeros(size(eta));
    ttr = zeros(size(eta));
    wtr = zeros(size(eta));

    eta = eta + abs(min(eta));

    [~,Plocs] = findpeaks(eta);
    [~,Tlocs] = findpeaks(-eta);
    dx = 10;
    ind = 1:length(eta);




    if numel(Plocs) == numel(Tlocs)
        h = eta(Plocs);
        t = eta(Tlocs);
        w = (diff(Plocs)*dx);
        lipks = ismember(ind,Plocs);
        lipks_w = ismember(ind,Plocs(1:end-1));
        htr(Plocs) = h;
        ttr(Plocs) = t;
        wtr(Plocs(1:end-1)) = w;

    elseif numel(Plocs) > numel(Tlocs)
        h = eta(Plocs(1:end-1));
        t = eta(Tlocs);
        w = (diff(Plocs(1:end-1))*dx);
        lipks = ismember(ind,Plocs(1:end-1));
        lipks_w = ismember(ind,Plocs(1:end-2));
        htr(Plocs(1:end-1)) = h;
        ttr(Plocs(1:end-1)) = t;
        wtr(Plocs(1:end-2)) = w;
    elseif numel(Plocs) < numel(Tlocs)
        h = eta(Plocs);
        t = eta(Tlocs(1:end-1));
        w = (diff(Plocs)*dx);
        lipks = ismember(ind,Plocs);
        lipks_w = ismember(ind,Plocs(1:end-1));
        htr(Plocs) = h;
        ttr(Plocs) = t;
        wtr(Plocs(1:end-1)) = w;
    else
        htr = NaN;
        ttr = NaN;
        wtr = NaN;
        lipks = NaN;
        lipks_w = NaN;
    end
    %calcualted by peaks

end

function [lipks,mal_m,mal_std,mal] = TwoDCrestFind_logical(eta)
    %function to find terminations in a logical matrix.

    for i = 1:size(eta,1)
        [~,locs]= findpeaks(eta(i,:));
        lipks(i,1:size(eta,2)) = ismember(1:size(eta,2),locs);
    end

    % for i = 1:size(eta,2)
    %     [~,locs]= findpeaks(eta(:,i));
    %     lipksW(1:size(eta,1),i) = ismember(1:size(eta,1),locs)';
    % end
    % lipks = lipksL | lipksW;
    lipks = bwmorph(lipks,'clean');
    lipks = bwmorph(lipks,'skel',inf);
    mal = regionprops(lipks,'MajorAxisLength','Centroid','PixelList');
    mal_m = mean([mal.MajorAxisLength]);
    mal_std = std([mal.MajorAxisLength]);
end

function [endpts] = TwoDCrestFind_logical_space(etas)
    %function to find terminations in a logical matrix.
    endpts = zeros(size(etas));

        for i = 1:size(etas,3)
            lipks = etas(:,:,i);
            lipks = lipks - repmat(smooth(mean(etas(:,:,1)),15)',size(etas,1),1);
            lipks = lipks>mean(lipks(:));
            lipks = bwmorph(lipks,'skel',inf);
            endpts(:,:,i) = bwmorph(lipks,'branchpoints') + bwmorph(lipks,'endpoints');
        end
end

function [endpts] = TwoDCrestFind_logical_time(etas)
    %function to find terminations in a logical matrix.
    for i = 1:size(etas,3)
        lipks = etas(:,:,i);
        lipks = lipks>mean(lipks(:));
        lipks = bwmorph(lipks,'skel',inf);
        endpts(i) = nnz(bwmorph(lipks,'branchpoints')) + nnz(bwmorph(lipks,'endpoints'));
    %     clts = regionprops(lipks-bwmorph(lipks,'branchpoints'),'MajorAxisLength');
    %     meancl(i) = mean([clts.MajorAxisLength]);
    end
end
  
 %% Utility functions
 
 function createLocalVar(inputStruct)
    %create local variables from the input structure.
    fns = fieldnames(inputStruct);

    for idx = 1:length(fns)
        assignin('caller',fns{idx},inputStruct.(fns{idx}));
    end
 end
 
 function [ina,da] = surf2theta(X,Y,Z,wd,method)
    % Find incidence and dip angle of surfaces
    %
    % if method == 1 : incidence angle is calculated per wind event direction,
    % wd (as long as wd is not empty)
    %
    % if method == 2 : dip azmuith calculation
    %
    % if run without wind direction, simply the surface (dip) azmuth is
    % returned in place of incidence angle.
    % function to calculate an array of incidence angles based on a single wind
    % vector and a 3d surface.



    if method == 1

        ina = zeros(size(Z));
        da = zeros(size(Z));

        if  ~isempty(wd)
            %this assumes that the user wants a incidence angle.

            wdx = cosd(wd);
            wdy = sind(wd);

            for i =1:size(Z,3)
                [nx,ny,nz] = surfnorm(X,Y,Z(:,:,i));
                hyp_ma = hypot(nx,ny);
                da(:,:,i) = 90-asind(nz./hyp_ma);  
                nx = nx./hyp_ma;
                ny = ny./hyp_ma;
                inadp= acosd(nx.*wdx(i) + wdy(i).*ny);
                %this is to create everything on the 0-90 domain.
                inadp(inadp>90) = 180-inadp(inadp>90);
                ina(:,:,i) = inadp;
            end

        elseif isempty(wd)
            %This assumes that the user wants a dip azmuith.


            for i =1:size(Z,3)
                [nx,ny,nz] = surfnorm(X,Y,Z(:,:,i));
                hyp_ma = hypot(nx,ny);
                da(:,:,i) = 90-atand(nz./hyp_ma); 
                nx = nx./hyp_ma;
                ny = ny./hyp_ma;
                ina(:,:,i) = rad2deg(cart2pol(-nx,-ny));
            end

        end

    elseif method == 2

          d=(-1:1)/dx; Zx=imfilter(Z,d,'circ'); Zy=imfilter(Z,d','circ');
          sz=size(Z); n=cat(3,-Zx,-Zy,ones(sz)); nmag=sqrt(sum(n.^2,3));
          n=bsxfun(@rdivide,n,nmag+~nmag);

    elseif (method ~= 1) && (method ~= 2)

        fprintf('Method choice is not valid\n\n');

    end
 end
 
 function [ina,da] = fast_surf2theta(X,Y,Z,wd,dx)

    % Find incidence and dip angle of surfaces!
    % if run without wind direction, simply the surface (dip) azmuth is
    % returned in place of incidence angle.
    % function to calculate an array of incidence angles based on a single wind
    % vector and a 3d surface.



    ina = zeros(size(Z));
    da = zeros(size(Z));


    dx = 10;

    if  ~isempty(wd)
    %this assumes that the user wants a incidence angle.

        wdx = cosd(wd);
        wdy = sind(wd);

        %assume equal grid

        if length(size(Z)) == 3
            [Fx,Fy,~]=gradient(Z,dx);
        else
            [Fx,Fy]=gradient(Z,dx);
        end

        N = Fx.^2 + Fy.^2 + 1;
        Nx = Fx./N;
        Ny = Fy./N;
        Nz = 1./N;
        ina = acosd(bsxfun(@times,cosd(wd),Nx) + bsxfun(@times,cosd(wd),Ny));
        da = asind(1-Nz);

    elseif isempty(wd)
        %This assumes that the user wants a dip azmuith.

        if length(size(Z)) == 3
            [Fx,Fy,~]=gradient(Z,dx);
        else
            [Fx,Fy]=gradient(Z,dx);
        end
            N = Fx.^2 + Fy.^2 + 1;
            Nx = Fx./N;
            Ny = Fy./N;
            Nz = 1./N;
            ina = rad2deg(cart2pol(Nx,Ny));
            da = asind(1-Nz);
    end
 end