%% A code for 2D Navier-stokes incompressible flow velocity/pressure formulation and staggered [with custom scenarios]
% created by Jamie Johns 2016 [optimized in 2018 but subject to improvement/optimization]

%sections of this code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ---> SECTION 1 - MASTER SECTION - SET PARAMETERS AND PERFORM MAIN CALCULATIONS FOR VELOCITY AND PRESSURE (Run this section first)
% ---> SECTION 2 - Visualise still frames of velocity field (used to create images in the examples of readme.pdf)
% ---> SECTION 3 - Record time evolution of velocity magnitude (no velocity field plot [streamslice])
% ---> SECTION 4 - Record time evolution of velocity magnitude AND velocity field [RECORDING TAKE LONGER THEN SECTION 3]
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% REQUIRED (DEPENDENT) MATLAB FILES; imagegrab.m, progressdisp.m, openvar.m, savevar.m,visualisation_still_frame.m
%                                     visualisation_record_video.m and visualisation_quick_animation.m


% Refer to "Readme.pdf" for example usage of this code as well as further information about the design of the code and useful sources for learning.


%This file (and previously other mentioned .m files) was can be found at;
%-->   https://github.com/JamieMJohns/Navier-stokes-2D-numerical-solve-incompressible-flow-with-custom-scenarios-MATLAB 



%% SECTION 1 - MASTER SECTION - SET PARAMETERS AND PERFORM MAIN CALCULATIONS
%TOTAL REVAMP OF GRID (to staggered half-griding)

%Initialize matlab for running code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nothing=input('<Press enter to begin calculations or cntrl+c to cancel>','s'); %This line is here to stop user from accidentally running this section twice (over-writting previous calculations)
    close all %close any open figures
    clear all %clear all variables in workspace
    clc %clear command window
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%Parameters for scenario (Modify these) ###############################################################################################################################################################################
    %set information about domain of simulation@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        SCENARIO='scenario_sphere6.png'; %<--- file (image) that is scenario for navier stokes fluid simulation
        domainX=2; % length of domain (x-axis) [unit: meters] 
        xinc=200; %number of nodes across x-component of domain (number of nodes from x=0 to x=domainX); where dx=domainX/xinc (=dy=dn)
        dt=0.0015; %set set delta time [unit: seconds]
        MI=3000; %number of time steps to perform calculations [time(at time step)=MI*dt]
        velyi=0; %y-component velocity of region with constant velocity (regions coloured red in scenario image)  [unit: meters/second]
                   %[velyi>0,velocity has vector -y with mag abs(velyi) and velyi<0, vel has vector +y with mag of abs(velyi)]
        velxi=1; %x-component velocity of region with constant velocity (regions coloured red in SCENARIO)   [unit: meters/second]
                   %[velxi>0,velocity has vector +x with mag abs(velxi) and velxi<0, vel has vector -x with mag of abs(velxi)]
        dens=1; %density  [unit: kg/m^3] , water(pure)=1000 blood=1025 air~1 
        mu=1/1000; %dynamic viscosity [kg/(m*s)]
        
    %Poisson Pressure solver parameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        error=0.001; %set tolerance of error for convergence poisson solver of Pressure field (good value to start with is 0.001; which is for most incompressible applications)
        MAXIT=1000; %maximum number of iterations allowed for poisson solver (increasing this will allow for further convergence of p solver)
        MINIT=1; %mininum number of iterations allowed for poisson solver (increasing this will allow for further convergence of p solver)
        % Note that: MINIT should be less than MAXIT
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    %save parameters $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    spacelim=5; %limit for hardrive space usage (in gigabytes) for externally saved data (x and y component velocity at each time step)
    ST=[100 100 500]; % FOR variables of dimensions of ST(1)xST(2) 
                      % save variable data for x and y component velocities
                      % in chuncks of files , each with ST(1)xST(2)xST(3)
                      % size matrix.........this reduces memory as only one file is openend at a time.
                      %(increasing ST(3) will reduce number of externally saved files [still using same amount of space)
                      %(decreasing ST(3) will reduce number of externally saved files [still using same amount of space)
                      %[Files; openvar.m and savevar.m are used]
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################################################################################################################################################################




% EVERYTHING BELOW HERE IS AUTOMATIC CALCULATIONS/OPERATIONS (NO EDITING REQUIRED)$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%Extract scenario from image@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    [velinit bounds outflow BLUE XI YI]=imagegrab(SCENARIO,xinc); %Uses imagegrab.m to extract scenario details from input image SCENARIO ["BLUE" is so far unused] 
    dn=domainX/XI; %dn=dx=dy (distance between nodes) [same result as domainY/YI]
    domainY=domainX.*(YI./XI); %domain length of y axis (domainY); 
                               % this is scaled with respect to;
                               % -> specified domainX
                               % -> and, w.r.t to dimensions of input image (Scenario)
    i=2:YI+1; %column index for calculations
    j=2:XI+1; %row index for calculations

    ts=round(prod(ST)./(XI.*YI)); %ts = number of timesteps to calculate u and v before saving as external file
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%initialize variables pressure and velocity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U=zeros(YI+2,XI+2); V=U; %initialise variable space for x and y velocity (staggered grid)
    P0=V;P01=V;P1=V;f=V; %initialise variable space for pressure (P0,P01,P1&f used at various stages in calculations)
    C=V; OF=V; %initialize variable space for variables that will help signal boundary conditions
    U(i,j)=(velinit).*velxi; % set velocity of constant velocity nodes for u
    v(i,j)=(velinit).*velyi; % set velocity of constant velocity nodes for v
    velx=0.5.*(U(2:YI+1,3:XI+2)+U(2:YI+1,2:XI+1)); % initialize matrix that will store unstaggered grid for u
    vely=0.5.*(V(3:YI+2,2:XI+1)+V(2:YI+1,2:XI+1)); % initialize matrix that will store unstaggered grid for u
    LRFX=velx; %initial velocity x [saved for,later, visualisation sections of code]
    LRFY=vely; %initial velocity y [saved for,later, visualisation sections of code]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Define variables that are used to signal particular boundary conditions
    % C(i,j)=1, if node i,j is not a solid wall node (=0 if otherwise)
        C(2:YI,2:XI)=1;C(i,j)=(bounds(i-1,j-1)==0);
        C(i,j)=(velinit(i-1,j-1)==0).*C(i,j)+(velinit==1);
    % CV(i,j)=1, when i,j is not solid wall,constant velocity region or outlfow node (=0 if otherwise)
        CV=C; %initial equal CV=C [CV=1 velocity time variant, CV=0 if velocity constant over time] 
        CV(i,j)=(velinit(i-1,j-1)==0); % set nodes of constant velocity C=0 (C=0=unchanged velocity)
        CV(i,j)=(bounds(i-1,j-1)==0).*CV(i,j); %if bound
    % CVC(i,j)=1, only if node i,j is in region of constant velocity  (=0 if otherwise)  
        CVC=V; %CVC=1 only solely for when a constant velocity condition is set.
        CVC(i,j)=velinit;
        CVC(i,j)=(velinit==1)+(velinit==0).*(outflow==0).*CVC(i,j);
    % OF(i,j)=1, if node is outflow node (=0 if otherwise)
        OF(i,j)=outflow;
%#####################################################################


% Redefine Outflow boundary nodes (if exist in scenario)$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 %This section is to reduce amount nodes which are outflow, keeping only exterior
 %nodes as outflow node (if defined): [only want exterior nodes to be outflow]
    %Set exterior nodes of OF to replicate closest interior node
        OF(:,end)=OF(:,end-1)==1; OF(:,1)=OF(:,2)==1;
        OF(end,:)=OF(end-1,:)==1; OF(1,:)=OF(2,:)==1;
    %Set any interior node,previously outflow node, to not be outflow node:
        BN=zeros(size(C)); BN(i,j)=1; %temporary variable interior nodes =1 , exterior equal zero
        OF(BN==1)=0; %set interior node OF to 0 (only exterior can be 1, if it is outflow node)
    %OF2(i,j)=1 if it has a neighbouring cell that is an outflow node (OF(i,j)=1):
        OF2=OF; %initialize OF2
        OF2(i,j)=(OF(i+1,j)+OF(i-1,j)+OF(i,j+1)+OF(i,j-1))~=0; %OF2(i,j)=1 if no OF=1 in neighbour cells
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%Define variable BC which determines orientation of Boundary with respect to node;
    BC=zeros(size(C)); %initialize boundary nodes
      %if BC(i,j)=2 , node (i,j+1) is either solid wall or constant velocity node [node (i,j) is "left" of solid wall or const vel node]  
        bc=(CV(i,j)).*(CV(i,j+1)==0);%if node i,j is a  free flow node but node i,j+1 is not (i.e- is a wall)
        BC(i,j)=(BC(i,j)==0).*((bc==1).*2+(bc==0).*BC(i,j))+(BC(i,j)~=0).*BC(i,j);
      %if BC(i,j)=4 , node (i,j-1) is either solid wall or constant velocity node [node (i,j) is "right" of solid wall or const vel node]       
        bc=(CV(i,j)).*(CV(i,j-1)==0); %if node i,j is a  free flow node but node i,j-1 is not 
        BC(i,j)=(BC(i,j)==0).*((bc==1).*4+(bc==0).*BC(i,j))+(BC(i,j)~=0).*BC(i,j);
      %if BC(i,j)=1 , node (i-1,j) is either solid wall or constant velocity node [node (i,j) is "above" solid wall or const vel node]     
        bc=(CV(i,j)).*(CV(i-1,j)==0); %if node i,j is a  free flow node but node i-1,j is not 
        BC(i,j)=(BC(i,j)==0).*((bc==1).*1+(bc==0).*BC(i,j))+(BC(i,j)~=0).*BC(i,j);
      %if BC(i,j)=3 , node (i+1,j) is either solid wall or constant velocity node [node (i,j) is "below" solid wall or const vel node]  
        bc=(CV(i,j)==1).*(CV(i+1,j)==0);%if node i,j is a  free flow node but node i+1,j is not
        BC(i,j)=(BC(i,j)==0).*((bc==1).*3+(bc==0).*BC(i,j))+(BC(i,j)~=0).*BC(i,j);
    BC(OF==1)=0; %exclude outflow nodes as being handled as solid wall/const vel boundary
    BC(OF2==1)=0; %exclude outflow nodes as being handled as solid wall/const vel boundary
    BC=(CV==1).*BC;
%#################################################################################################


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%determine divisor for pressure poisson equation$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    CP=zeros(size(P0)); %initial matrix to store divisor coefficients in poisson solver
    CP(i,j)=(C(i+1,j)==1)+(C(i-1,j)==1)+(C(i,j-1)==1)+(C(i,j+1)==1); %only solid walls are zeros  %<<<<<<<<<<<<<<<<<<<<<<<<<typically using
    %-->CP(i,j)=4 if node i,j not next to solid wall (interior node)
    %-->CP(i,j)=3 if node i,j next solid wall (above/down/left/right)
    %-->CP(i,j)=2 if node i,j next corner of solid wall (corner node)
    %-->node i,j is surrounded by solid wall nodes below will be zero (but correct to 1)
    CP(CP<=1)=4; % if CP(i,j)==0, set to 4 to avoid division by zero in poisson equation (happens if node i,j is "inside" solid wall)
    CP(OF2==1)=4; % if node i,j is neighbouring outflow node [OF2(i,j)=1]; set divisor to 4
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


% DETERMINE MEMORY USAGE$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    if exist([pwd '\temporary_NS_velocity'],'dir')==0 %if directory does not exist (external save of u and v)
        mkdir([pwd '\temporary_NS_velocity']) %create directory to store u and v
    end
    save([pwd '\temporary_NS_velocity\' '00_TEMPORARY_FILES_4_NAV_STK_CODE' '.mat' ],'XI') %save random file to directory  (testing)
    % 1 kilobyte = 2^10 bytes , 1 megabyte = 2^10 kilobytes , 1 gigabyte = 2^10 megabytes
    % ans=1 (variable) , ans=8 bytes  ; ans=ones(2,2,3) ,ans=2*2*3 bytes
    %matrix of size a x b x c = a*b*c*8 bytes = a*b*c*8./(C^3) gigabytes ; C=2^10
    % FILES THAT WILL BE SAVE velx and vely; so if mesh = 200 x 200 and MI=1000
    % total saved on hardrive space = 8*2*200*200*MI/(C^3)
    dirr=[pwd '\temporary_NS_velocity\']; %record directory path for print (below)
    spaceused=1.0889.*(8*XI*YI*MI/(10^9)); %determine approximate space used
    fprintf('\n It is estimated that approximately %.4f gigabytes of hard drive space \n',spaceused);
    fprintf('These temporary files will be saved in: \n');
    fprintf('%s \n',dirr);
    fprintf('For each velocity (x and y) there will be %.0f files \n',ceil(MI./ts));
    fprintf('therefor %.0f files with average size of %.6f gigabytes \n',2.*ceil(MI./ts),1.0889.*(8*XI*YI*ts)/(2.*10^9));
    if spaceused>spacelim; %if estimated space used is more than limitations set by user (flash warning)
       fprintf('\nWarning!!!!\n')
       fprintf('The above estimated hard drive space usage (%.4fgigabytes)\n',spaceused);
       fprintf('is more than the limit set by user of %.4fgigabytes.\n',spacelim);
       userpromp=input('<Press enter to continue or cntrl+c to cancel calculations>','s');
       commandwindow %bring command window to foreground
    end
    ts=round(prod(ST)./(XI.*YI)); %parameter that determines number of timesteps of calculating
                                  % u and v before they are externally saved to a file (on hardrive)
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



%Anonymous functions ##########################################################################################################
    dxuu=@(a,b,i,j) (((a(i,j+1)+a(i,j)).*0.5).^2-((a(i,j)+a(i,j-1)).*0.5).^2)./dn; %finite difference derivative of U*U w.r.t x  (staggered grid) 
    dyuv=@(a,b,i,j) ((a(i+1,j)+a(i,j)).*(b(i,j+1)+b(i,j)).*0.25-(a(i,j)+a(i-1,j)).*(b(i-1,j+1)+b(i-1,j)).*0.25)./dn; %finite difference derivative of U*V w.r.t y  (staggered grid)
    dyvv=@(a,b,i,j) (((a(i+1,j)+a(i,j)).*0.5).^2-((a(i,j)+a(i-1,j)).*0.5).^2)./dn; %finite difference derivative of V*V w.r.t y  (staggered grid)
    dxuv=@(a,b,i,j) ((a(i+1,j)+a(i,j)).*(b(i,j+1)+b(i,j)).*0.25-(a(i,j-1)+a(i+1,j-1)).*(b(i,j)+b(i,j-1)).*0.25)./dn; %finite difference derivative of U*V w.r.t x  (staggered grid) 
    DX2=@(a,i,j) (a(i-1,j)+a(i+1,j)+a(i,j+1)+a(i,j-1)-4.*a(i,j))./(dn.^2); %finite difference for laplace of operator in two dimensions (second deriv x + second deriv y;staggered grid)
    dx0=@(a,i,j) (a(i,j+1)-a(i,j))./dn; % first order derivative finite difference used for pressure w.r.t x (un-staggered grid)
    dy0=@(a,i,j) (a(i+1,j)-a(i,j))./dn; %first order derivative  finite difference used for pressure w.r.t x (un-staggered grid)
%##############################################################################################################################



% Ininitialise parameters and variables used in calculations$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    [x,y]=meshgrid(linspace(0,domainX,XI),linspace(0,domainY,YI)); % x and y coordinates [used in visualisation sections]
    mvel=max(max(sqrt(velx.^2+vely.^2))); %initial maximum velocity (used for plotting later in sections 3 and 4)
    i=2:YI+1; %column index for calculations
    j=2:XI+1; %row index for calculations
    L=mu/dens;  %kinematic viscosity (
    TIME=0; % paramer which keeps track of simulated time
    docalc=1; %parameter that determines if calculations should continue (1=yes) 
    T=0;  %parameter which keeps track of total number of time steps calculated
    T2=0; %parameter that determines when set of u and v calculations should be save externally
    TL=0; %parameter that will record simulation time of last file saved
    vt=0; %parameter for function which prints progress of calculations in 10% increments
    nsave=0; %parameter that keeps track of file number for external saving of u and v
    slcv=0; %parameter for slow convergence message
    PIT=0; %array that will keep record of number of iterations during each time pressure poisson equation is solved
    tic; %start timer for calculations
    TLS=1;%time step at last save frame;
    T0=1; %time step of last calculate frame
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%CP=4.*ones(size(CP));
% cv2 from last happy!!!!!!!!!!!!!!
CA=zeros(size(CV));

CA(i,j)=(CVC(i,j)==1).*((CV(i+1,j)+CV(i-1,j)+CV(i,j+1)+CV(i,j-1))~=0);
CA=CA+CV;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% THE ACTUAL CALCULATIONS (Calculating TIME-DEPENDENT VELOCITY AND PRESSURE)@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

PA=P0;%temp use

while docalc==1; %WHILE CALCULATIONS ARE BEING RUN
    vt=progressdisp(T,MI,'Calculations','time-steps',vt,toc);  %function prints progress of calculations (string in command window) (by increment of 10%) [uses function progressdisp()]


% Specify velocities on constant velocity region @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    % (\[ CVC(i,j)=1 if node i,j is in constant velocity region]
      U(i,j)=(~CVC(i,j)).*U(i,j)+(CVC(i,j)).*velxi; %set x-velocity to velxi if node is in contant vel B.C  
      V(i,j)=(~CVC(i,j)).*V(i,j)+(CVC(i,j)).*velyi; %set y-velocity to velyi if node is in contant vel B.C      
 %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    

if T~=0 %if not instant of calculation
% Apply boundary conditions for nodes next to solid wall or next to constant velocity $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        % Non-slip boundary condition; e.g - B(i,j)=1 & CVC(i-1,j)=0  [wall exist above i,j and is non-moving solid; parallel vel is zero and tangent is reflect vel from i+1,j]
        % free-slip boundary condition; e.g - B(i,j)=2 & CVC(i,j+1)=1  [wall exist right of i,j and is moving wall (interpolate velocity at i,j); parallel is mean of vel for nodes i,j+1 and i,j-1]
            %If node node i,j is "below" wall or const vel [BC(i,j)=1]####################
             U(i,j)=(BC(i,j)~=1).*U(i,j)+(BC(i,j)==1).*((CVC(i-1,j)).*0.5.*(U(i-1,j)+U(i+1,j))); %velocity parallel with wall 
             V(i,j)=(BC(i,j)~=1).*V(i,j)+(BC(i,j)==1).*((CVC(i-1,j)).*0.5.*(V(i-1,j)+V(i+1,j))+(~CVC(i-1,j)).*(-V(i+1,j)).*2/3); %velocity tangent to wall
            %If node node i,j is "above" wall or const vel [BC(i,j)=3]####################
             U(i,j)=(BC(i,j)~=3).*U(i,j)+(BC(i,j)==3).*((CVC(i+1,j)).*0.5.*(U(i-1,j)+U(i+1,j))); %velocity parallel with wall (zero if wall is non-moving solid)
             V(i,j)=(BC(i,j)~=3).*V(i,j)+(BC(i,j)==3).*((CVC(i+1,j)).*0.5.*(V(i+1,j)+V(i-1,j))+(~CVC(i+1,j)).*(-V(i-1,j)).*2/3); %velocity parallel to wall Boundary (zero if wall)                
            %If node node i,j is "left" of wall or const vel [BC(i,j)=2]####################
             V(i,j)=(BC(i,j)~=2).*V(i,j)+(BC(i,j)==2).*((CVC(i,j+1)).*0.5.*(V(i,j-1)+V(i,j+1))); %velocity parallel with wall 
             U(i,j)=(BC(i,j)~=2).*U(i,j)+(BC(i,j)==2).*((CVC(i,j+1)).*0.5.*(U(i,j-1)+U(i,j+1))+(~CVC(i,j+1)).*(-U(i,j-1)).*2/3); %velocity tangent to wall                     
            %If node node i,j is "right" of wall or const vel [BC(i,j)=4]####################
             V(i,j)=(BC(i,j)~=4).*V(i,j)+(BC(i,j)==4).*((CVC(i,j-1)).*0.5.*(V(i,j-1)+V(i,j+1))); %velocity parallel with wall 
             U(i,j)=(BC(i,j)~=4).*U(i,j)+(BC(i,j)==4).*((CVC(i,j-1)).*0.5.*(U(i,j-1)+U(i,j+1))+(~CVC(i,j-1)).*(-U(i,j+1)).*2/3); %velocity tangent to wall  
   
 end   

    %apply neumann boundary condition for outflow nodes @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    % if node OF(i,j)=1 , where i,j is exterior node, u(i,j) and v(i,j) are equal to velocity at nearest nodes (which are not outlow nodes)
     % e.g- if OF(1,2:end-1)=1, then set; U(1,2:end-1)=U(2,2:end-1) and V(1,2:end-1)=V(2,2:end-1)
        % U velocity (handle exterior nodes)
            U(:,end)=(OF(:,end)).*U(:,end-1)+(~OF(:,end)).*U(:,end);%right boundary
            U(:,1)=(OF(:,1)).*U(:,2)+(~OF(:,1)).*U(:,1); %left boundary
            U(1,:)=(OF(1,:)).*U(2,:)+(~OF(1,:)).*U(1,:); %Top boundary
            U(end,:)=(OF(end,:)).*U(end-1,:)+(~OF(end,:)).*U(end,:);%Bottom boundary
        % V velocity (handle exterior nodes)
            V(:,end)=(OF(:,end)).*V(:,end-1)+(~OF(:,end)).*V(:,end);%right boundary
            V(:,1)=(OF(:,1)).*V(:,2)+(~OF(:,1)).*V(:,1);%left boundary
            V(1,:)=(OF(1,:)).*V(2,:)+(~OF(1,:)).*V(1,:);%Top boundary
            V(end,:)=(OF(end,:)).*V(end-1,:)+(~OF(end,:)).*V(end,:);%Bottom boundary
        % V velocity (handle exterior (corner) nodes)
            V(1,1)=(OF(1,1)).*V(2,2)+(~OF(1,1)).*V(1,1); %Top left corner
            V(1,end)=(OF(1,end)).*V(2,end-1)+(~OF(1,end)).*V(1,end); %Top right corner
            V(end,1)=(OF(end,1)).*V(end-1,2)+(~OF(end,1)).*V(end,1); %Bottom left corner
            V(end,end)=(OF(end,end)).*V(end-1,end-1)+(~OF(end,end)).*V(end,end);%Bottom right corner
        % U velocity (handle exterior (corner) nodes)
            U(1,1)=(OF(1,1)).*U(2,2)+(~OF(1,1)).*U(1,1); %Top left corner
            U(1,end)=(OF(1,end)).*U(2,end-1)+(~OF(1,end)).*U(1,end); %Top right corner
            U(end,1)=(OF(end,1)).*U(end-1,2)+(~OF(end,1)).*U(end,1); %Bottom left corner
            U(end,end)=(OF(end,end)).*U(end-1,end-1)+(~OF(end,end)).*U(end,end);%Bottom right corner
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


uc=U; %record current x-vel (for use later on to enforce constant x vel conditions)
vc=V; %record current x-vel (for use later on to enforce constant y vel conditions)


% Solve temporary velocities#########################################################################################################
        %calculate temporary velocity using convective and diffusive terms of momentum equations
        if T~=0 %if not first instant of calculation
            %if node is not solid wall or constant velocity region (CV=1) , update velocity
            %and, also, if node is not next to boundary(i.e - BC==0), update velocity [nodes i,j for B(i,j)~=0 already had velocities determined when apply B.C (above)]
            U(i,j)=U(i,j)+(BC(i,j)==0).*(CV(i,j)==1).*dt.*(L.*(DX2(uc,i,j))-(dxuu(uc,uc,i,j)+dyuv(uc,vc,i,j)));%update x-velocity 
            V(i,j)=V(i,j)+(BC(i,j)==0).*(CV(i,j)==1).*dt.*(L.*(DX2(vc,i,j))-(dxuv(uc,vc,i,j)+dyvv(vc,vc,i,j)));%update y-velocity          
        end
%####################################################################################################################################

      
%Finite difference solve for poisson pressure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %continuity equation is used in RHS side of poisson equation (f) to enforce incompressible flow
    errach=(T==0); %parameter that determine pressure poisson solution has converged
                       %or surpassed maximum allowed iterations (MAXIT)
                       %also, calculation not calculated for first instance (errach=1 if T=0)
    mite=0; %counts number of iterations run for poisson equation
    cp=(T~=1).*1+(T==1).*10; %multiply maximum allowed iterations, for poisson equation solve, by 10 to allow for smoothing at first step of calculation (T=1)
            f(i,j)=(OF2(i,j)==0).*(CA(i,j)==1).*(dn/dt).*(V(i,j)+U(i,j)-V(i-1,j)-(U(i,j-1))); % [right hand side of poisson equation] 
            while errach==0 %while solution has not converged for poisson equation
               P1(i,j)=0.75*(P0(i+1,j)+P0(i-1,j)+P0(i,j-1)+P0(i,j+1)-f(i,j))./CP(i,j)+0.25.*P0(i,j); %Finite difference solve of poisson pressure equation (for pressure field)   
              %  P1(i,j)=0.01*(P0(i+1,j)+P0(i-1,j)+P0(i,j-1)+P0(i,j+1)-f(i,j))./CP(i,j)+0.99.*P0(i,j); %Finite difference solve of poisson pressure equation (for pressure field)  
                P1(OF==1)=0; %if node is an outflow condition; set pressure to zero                    %using relaxtion method  A=0.75 and 1-A=0.25 (last term)
                P1(CA==0)=0; %if node is solid wall node (still or moving wall); set pressure to zero
                errach=(max(max(abs(P1-P0))))<=error; %determine if max difference between next and current pressure
                                                     %field satifies allowed tolerance for error
                                                     %[If yes, errach= 1 = solution converged (exit while loop)]
                mite=mite+1; %add to count of iterations of poisson equation solution
                PFF=abs(P1-P0);
                P0=P1; %update current pressure field with next step
               % PA(:,:,size(PA,3)+1)=P0;%temp use
                if mite>=(MAXIT*cp); %if iterations of poisson solver surpases max number of allow iterations
                 errach=1; %set errach =1 (exit while loop) regardless of whether tolerance of error has not been reached
                end
                if mite<MINIT %if iterations of poisson solver is less than mininum required number of iteration (MINIT)
                     errach=0; %set errach =0 (don't exit while loop) regardless of whether tolerance of error has been reached
                end

                %Prompt to indicate to user that poisson solver might be taking long time to solve######################################
                        if T>100 &&  slcv==0 && (mean(PIT)>(MAXIT)/2); %if mean number iterations after 100 frames is greater than half MAXIT
                            fprintf('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                            fprintf('Slow convergence for Pressure poisson solver (P.P) was detected!\n')
                            fprintf('To increase speed of convergence you could stop calculations and modify some\n')
                            fprintf('parameters such as;\n')
                            fprintf('-> decrease maximum allowed iterations P.P ("MAXIT")\n')
                            fprintf('-> tolerance of error for pressure poisson ("error")\n')
                            fprintf('-> delta time (dt)\n')
                            fprintf('-> resolution of scenario ("MI")\n')
                            fprintf('-> reynolds number (RE)\n')
                            fprintf('-> velocity magnitude at boundary conditions ("velxi") and/or ("velyi")\n')
                            fprintf('In the mean,time calculations will run until manually stopped\n')
                            fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                            slcv=1; %update prompt parameter such that prompt will not come up a second time
                        end
                %########################################################################################################################

            end
            PIT(size(PIT,2)+1)=mite; %record number of iterations ran for poisson solver (DEBUG USE)
       %Prompt to warn user that pressure solution has diverged (and also forcibly stop calculations [since no good will come from further running code])
        if sum(sum(isnan(P0)))~=0 || sum(sum(isinf(P0)))~=0 % if Pressure has infinite or nan value
          fprintf('\n Error! The solution for pressure has diverged (going to inf or nan)\n')
          fprintf('Try modifying dt,mu,xinc,domainX,velxi/velyi and try again!\n')
          force_exit_code=intentional_error; % this line intentionally causes error to stop code (remove this line if code won't run at all)
        end
       %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CORRECT Velocities BY ADDING PRESURE GRADIENT$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
         if T~=0 %if not first instant being recorded (where time=0 seconds)
            %As with calculating temporary velocity (a previous section);
            %   If node is not solid wall or constant velocity region (CV=1) , update velocity
            %   And, also, if node is not next to boundary(i.e - BC==0), update velocity [nodes i,j for B(i,j)~=0 already had velocities determined when apply B.C (above)]
            U(i,j)=U(i,j)+(OF(i,j)==0).*(BC(i,j)==0).*(CV(i,j)==1).*dt.*(-(dx0(P0,i,j)./dens)); %update x-velocity
            V(i,j)=V(i,j)+(OF(i,j)==0).*(BC(i,j)==0).*(CV(i,j)==1).*dt.*(-(dy0(P0,i,j)./dens)); %update y-velocity 
         end
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
         
 
%if T==1;toc; falhg=adgjlkadgh;end;
% Obtaining velocities on unstaggered grid and saving temporary files###################################################################################################################################

    % obtain velocities for unstaggered grid @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
            u1=0.5.*(U(2:YI+1,3:XI+2)+U(2:YI+1,2:XI+1)); %x-comp velocity of unstaggerd grid (velx): velx(i,j)=(U(i,j+1)+U(i,j))/2 [U=staggered grid x vel]
            v1=0.5.*(V(3:YI+2,2:XI+1)+V(2:YI+1,2:XI+1)); %y-comp velocity of unstaggerd grid (vely): vely(i,j)=(V(i+1,j)+V(i,j))/2 [V=staggered grid y vel]
            velx(i-1,j-1,T2+1)=(CV(i,j)).*u1+(~CV(i,j)).*uc(i,j); % record unstaggered velocity to list velx
            vely(i-1,j-1,T2+1)=(CV(i,j)).*v1+(~CV(i,j)).*vc(i,j); % record unstaggered velocity to list vely
            TIME=TIME+dt; %update time (simulation time)
            T2=T2+1; %update parameter which determines when to externally save set calculations (x,y velocity)
            T=T+1; %update parameter which is total number of time steps
            T0=T; %time of last calculated frame
            LRFX=velx(:,:,T2); %record last calculated instance of velocity x [for use in visualization if calculations are halted before completeion]
            LRFY=vely(:,:,T2); %record last calculated instance of velocity y [for use in visualization if calculations are halted before completeion]
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    % determine if current set of x and y velocity should be saved$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if T2>= ts %=True if current x,y velocity needs to be saved
            nsave=nsave+1; %update file save number
            mvel1=max(max(max(sqrt(velx.^2+vely.^2)))); %current maximum recorded velocity
            mvel(mvel1>=mvel)=mvel1; %if current max vel is greater than max recorded vel, then update
            storevar(velx,'NSTOKES_TEMP_vx_',nsave) %save x velocity as file number nsave
            storevar(vely,'NSTOKES_TEMP_vy_',nsave) %save y velocity as file number nsave
            TL=TIME; %record last simulation time
            TLS=T; %record timestep of last saved frame
            T2=0;  %reset save parameter to zero   
        end
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


    %determine if calcuations are finished%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if T>=MI %if maximum allowed number of timesteps is reached
            docalc=0; %set parameter to finish calculations (while loop does not continue)
            if T2~=0 %if last set of calculations was not already saved
                %then save final set of calculations (keeping only vely and velx calculated in last set (1 to T2);
                nsave=nsave+1; %update file save number
                storevar(velx(:,:,1:T2),'NSTOKES_TEMP_vx_',nsave)  %save x velocity as file number nsave
                storevar(vely(:,:,1:T2),'NSTOKES_TEMP_vy_',nsave)  %save y velocity as file number nsave
                TL=TIME; %record last simulation time
                TLS=T; %record timestep of last saved frame
            end
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
%#############################################################################################################################################################################################################
    
end %if docalc=1, restart while loop [continue calculations
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

fprintf('\nCalculations are 100%% complete! at %.4f seconds [%.0f of %.0f time-steps]\n',toc,T,MI)

%% Section 2 - View velocity field at last step
%In this section (after running first section), three figures will be created;
%   figure(1); visualisation of velocity at last calculated instant
%   figure(2); visualisation of velocity at last instant that was saved as an external file
%   figure(3); visualisation of velocity at a particular instant specified by user (which is,below; paramter "TF")
%Note: if calculations were not haulted before complete (in section 1), figure(1) and figure(2) in this section will be the same as each other

%Also note that for figures (in this section), if the figure title appears to 
%be cut off, try stretching the figure [this should show all text]
%Options for plotting@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%[note: the parameters,below, are used in "visualisation_still_frame.m" (which is also execute below)]
    velfield=1; %1=show velocity field (indicated by streamslice() function 
    arrowdensity=3; %density of arrows for streamslice()                    
    arrowwidth=1.50; %width of arrows for streamslice()
    TF=1300; %specific instant to plot for figure(3) ("Target Frame") [should be equal to, or less than, MI (section 1 parameter)]/
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


run visualisation_still_frame.m   %execute script (below) which will produce the three previously mentioned figures

%% SECTION 3- quick visualisation of velcoity (No video recording option)
%This section is intended for quick (animated) visualisation of the calculated
%velocities (showing how the velocity changes w.r.t time).

%Note: The real time speed of the animated figures will depend on the user's computer specifications (i.e - RAM,GPU,CPU, etc...)
%      as well as depend on some settings (parameters) that can be specified below. 

%parameters$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%[note: the parameters,below, are used in "visualisation_quick_animation.m" (which is also execute below)]
    velfield=0; %if =1, indicate velocity direction with streamslice() function [note that this could slow down rate of animation speed]
    arrowdensity=2; %density of arrows for streamslice() [default = 2]
    arrowwidth=1; %width of arrows for streamslice() [default = 1]
    title_detailed=0; % if =1, title in figure has time-dependent details about scenario [note that this could slow down rate of animation speed]
    colourbar=1; %if =1, in figure; show colourbar (for velocity magnitude) [note that this could slow down rate of animation speed]
    capvel=0; %if =1, colour bar ranges between 0 and maximum recorded vel mag (section 1), else; colourbar axis limits can stretch and shrink if new velocity is outside limit of colourbar
    tinc=1; %number of frames to skip plotting (i.e - if tinc=20 , frames to be plotted are 1,21,41,...etc)
              %[increasing tinc will increase animation speed but may miss out finer details)
    fntsize_title=11; %font size of title in figure
    fntsize_axis_label=13; %font size of axis label in figure
    fntsize_axis_label_colorbar=12; %font size of label on colour bar
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

run visualisation_quick_animation.m  %execute script (below) which will produce the three previously mentioned figures
%% title image
xlabel('X(meters)','fontsize',45)
ylabel('Y(meters)','fontsize',45)
ylabel(cb,'Velocity Magnitude (m/s)','fontsize',40) 
title('','fontsize',20)
%% SECTION 4 - visualisation of velocity including video recording

%This section is the same as section section 3, however is for recording videos (in .avi format from it).
%When you run this section, you will be asked to stretch the video to a desired resolution for the 
%recorded video [this will be actual resolution of the video (smaller resolution=faster video recording = less detail)].

%And, as a video is being recorded from the animated figure; the figure should be kept in the foreground
%of the window during the recording. (anything in front of the figure will be recorded into the video).

%Also note that the real-time speed of the animated figures (and time taken to record the video)
%will depend on the user's computer specifications (i.e - RAM,GPU,CPU, etc...) as well as depend on some settings 
%(parameters) that can be specified below.
 


%Parameters$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    %Figure (animation) parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[note: the parameters,below, are used in "visualisation_quick_animation.m" (which is also execute below)]
        velfield=1; %if =1, indicate velocity field with streamslice() [note that this could slow down rate of animation speed]
        arrowdensity=4; %density of arrows for streamslice() (if velfield=0, this parameter does not apply) [default: 2]
        arrowwidth=1; %width of arrows for streamslice() (if velfield=0, this parameter does not apply) [default: 1]
        colourbar=1; %if =1, in figure; show colourbar (for velocity magnitude) [note that this could slow down rate of animation speed]
        capvel=1; %if =1, colour bar ranges between 0 and maximum recorded vel mag (section 1), else; colour bar updates on each frame (depending on max/min velocity magnitudes)
        title_detailed=1; % if =1, title in figure is has details about scenario [note that this could slow down rate of animation speed, ALTHOUGH, turning this on will help user keep track of progress of video recording]
        fntsize_title=11; %font size of title in figure
        fntsize_axis_label=13; %font size of axis label in figure
        fntsize_axis_label_colorbar=12; %font size of label on colour bar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Video recording parameters###########################################
        RECORD=1; %if =1, record animated output to video file (.avi)
        videoname='Video5'; %filename of video being recorded (i.e - if videoname='myvideo'; filename is "myvideo.avi")
        timetarg=15; % desired length of output video in seconds
        framespersec=30; %frames per second for output video
        %Note: depending on the of number calculated time-steps of velocity (section 1)
        %      and required number of frames for the recorded video [total needed frames =timetarg*framespersec]
        %      some time-steps may be skipped, or recorded multiple times, to satisfy desired timelength (timetarg) and frame-rate (framespersec).
    %######################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

run visualisation_record_video.m  %execute script (below) which will produce the three previously mentioned figures

