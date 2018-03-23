%created by Jamie Johns 2018
%NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%This .m file is accessed and executed within "main.m" ( specifically"section 3" after running
%"section 1").
%And running this .m file, independently, may not work without first running the previously mentioned "sections" of main.m
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

[x,y]=meshgrid(linspace(0,domainX,XI),linspace(0,domainY,YI)); %  x and y coordinates for plotting
clc %clear command window
close all%close any open figures

% Plot figure at initial state to initialise plot objects##################################
figure(1) 
    hold on %hold subsequently plotted objects without deleting previous object
    velx=openvar('velx','NSTOKES_TEMP_vx_',1); %get u from first externally saved file
    vely=openvar('vely','NSTOKES_TEMP_vy_',1); %get v from first externally saved file
    vel_mag=pcolor(x,y,sqrt(velx(:,:,1).^2+vely(:,:,1).^2)); %plot pcolor object to display velocity magnitude
    shading interp %smooth shading
    WALLS=-10.*ones(size(x)); %matrix that will represtent solid walls
    WALLS(bounds==0)=nan; %if matrix element is not solid wall node (set to nan so is invisible)
    SOLIDWALL=surf(x,y,WALLS); %create surf object that will display solid wall nodes (coloured black)
    set(SOLIDWALL,'facecolor',[0,0,0]) %set surf object to black
    xlabel('X(meters)','fontsize',fntsize_axis_label) %set label for x axis
    ylabel('Y(meters)','fontsize',fntsize_axis_label) %set label for y axis
    view(0,-90) %flip plot view (vertically)
    xlim([0 domainX]) %set x-axis plot limits
    ylim([0 domainY]) %set y-axis plot limits
    axis image %set image axis to stretch to image specifications
     if colourbar==1 %if colour bar to be displayed
        cba=colorbar; %create colour bar
        ylabel(cba,'Velocity magnitude (m/s)','fontsize',fntsize_axis_label_colorbar) %add y label to colour bar
      if capvel==1; caxis([0 mvel]); end; %if colour bar is kept constant w.r.t maximum calculated velocity (over all time steps)
     end
%################################################################################################
 cax=caxis;
%determine list of time-steps to plot u,v@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
tinc(tinc<1)=1; %mininum increement of time step is one
if tinc==1 % if time steps to plot from 1 to final frame is in steps of 1
    tf=1:MI;
else %else, if in steps of "tinc"
    tf=[1:tinc:MI MI];   %list of frames from 1 to MI in steps of "tinc"
end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%Below is part of string used in title in figure (STATIC PART OF STRING)
strmg1=sprintf('Rendered in Matlab with code created by Jamie M. Johns\n density=%.2fkg/m^3; \\mu=%.4fkg/(m*s); dt=%.4fs, resolution:%.0fx%.0f [for calculations]',dens,mu,dt,XI,YI);
title(sprintf('%s',strmg1),'fontsize',fntsize_title) %initial title for animation

T2=1; %parameter which keeps track when next externally saved data (for u and v) should be loaded
ns=1; %parameter which keeps track which externally saved data to load
vel_field=[]; %initialise parameter that will be object for plotted streamslice (velocity field)
M=sqrt(velx.^2+vely.^2); %velocity magnitude
fcnt=1;%frame counter
caxis([0 mvel1])
vt=0; %parameter for function which prints progress of calculations in 10% increments
tic
for T=tf %for each time step to be plotted
    vt=progressdisp(fcnt,length(tf),'plotting','frames',vt,toc);  %function prints progress of calculations (by increment of 10%)
    set(vel_mag,'cdata',M(:,:,T2)); %update colour data to current velocity magnitudes

    if velfield==1; %if indicate velocityfield (using stream slice)$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if ~isempty(vel_field); delete(vel_field); clear hd; end; %delete vel_field (from last time-step)
        vel_field=streamslice(x,y,velx(:,:,T2),vely(:,:,T2),arrowdensity); %streamslice and set arrow/line density
       set(vel_field,'color',[0.0 0.0 0.0],'linewidth',arrowwidth) %set colour (to black) and width of lines
    end%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    %update title to current time step information:c
    if title_detailed==1; %if using detailed title
        title(sprintf('%s\n time=%.4f seconds; frame:%.0f/%.0f',strmg1,T.*dt,fcnt,length(tf)),'fontsize',fntsize_title)
    end
    drawnow %update plot objects visually on figure
    fcnt=fcnt+1; %update frame counter
    T2=T2+tinc; %update parameter which keeps track of time step w.r.t to current loaded file
    if T2>=ts %determine if next file needs to be loaded
        ns=ns+1;  %file number to load (next set of calculations
        if ns<=nsave %if not at maximum number of saved files
           velx=openvar('velx','NSTOKES_TEMP_vx_',ns); %load next set of saved u
           vely=openvar('vely','NSTOKES_TEMP_vy_',ns); %load next set of saved v
           M=sqrt(velx.^2+vely.^2); %calculate velocity magnitude (for next set of u an v)
        end
       T2=1+(T2-ts); %update time step parameter w.r.t to newly loaded files
    end
end 
fprintf('\nCalculations are 100%% complete! at %.4f seconds\n',toc)