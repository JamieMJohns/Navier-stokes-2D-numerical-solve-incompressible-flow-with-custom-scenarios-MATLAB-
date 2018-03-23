%created by Jamie Johns 2018
%NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%This .m file is accessed and executed within "main.m" ( specifically "section 2" after running
%"section 1").
%And running this .m file, independently, may not work; especially without first running the previously mentioned "sections" of main.m
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%below is static part of string for titles in figure
strmg1=sprintf('Rendered in Matlab with code created by Jamie M. Johns\n density=%.2fkg/m^3; \\mu=%.4fkg/(m*s); dt=%.4fs, resolution:%.0fx%.0f [for calculations]',dens,mu,dt,XI,YI);
close all %close all open figures
clc %clear command window
[x,y]=meshgrid(linspace(0,domainX,XI),linspace(0,domainY,YI)); % x and y coordinates for u and v



% PLOT VELOCITY OF LAST CALCULATED FRAME$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

figure(1) %plot of velocity for last calculated time step (of section 1)
    %from section 1: LRFX= last calculated u and LRFX= last calculated v
    LRFX((bounds)==1)=nan; % set nodes inside solid wall region to nan (does not show up in velocity plot)
    LRFY((bounds)==1)=nan; % set nodes inside solid wall region to nan (does not show up in velocity plot)
    F=flipud(sqrt(LRFX.^2+LRFY.^2)); % velocity magnitude
    hold on %hold all objects created on figure (create new object without deleting previous)
    pcolor(x,y,F); %plot velocity magnitude for each node w.r.t x and y position
    xlabel(sprintf('x(meters)\n[time-step:%.0f of %.0f (last calculated timestep)]',T0,MI)) %label x axis
    ylabel('y(meters)') % label y axis
    title(sprintf('%s\nsimulation time=%.4fseconds',strmg1,T0.*dt)) %title for figure (show time to 4 decimal places [%.4f])
    cb=colorbar; %create colorbar object with reference "cb"
    cba=caxis; %record current min/max values of colour bar (and hence min/max magnitude of velocity)
    ylabel(cb,'Velocity (m/s)') %add ylabel to colour bar
    shading interp %use smooth shading for pcolor (comment out this line to see the difference).
    if velfield==1 %if velocityfield is to be indicated (by streamslice() function)
        hd=streamslice(x,y,flipud(LRFX),-flipud(LRFY),arrowdensity); %plot streamslice (indicate velocity field)
        set(hd,'color',[0.0 0.0 0.0],'linewidth',arrowwidth) %set colour and line width of arrows in streamslice
    end
    sld=surf(x,y,ones(size(x)).*-1); %create foreground object (used to highlight solid boundaries with black colour)
    shading interp %reinforce smooth shading for pcolor and surf objects
    axis image %stretch axis to scale image
    set(sld,'facecolor',[0.0 0.0 0.0]) %set colour of sld to black
    caxis([cba(1) cba(2)]) %reset colorbar axis to only consider F (and not sld)

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


%PLOT VELOCITY OF LAST EXTERNAL SAVED FILE (AND LAST TIME STEP) ####################################################################################
if nsave~=0; %if external files were saved (nsave will not equal zero)
    figure(2) %plot of velocity for the last timestep saved in external file
        hold on %hold all objects created on figure (create new object without deleting previous)
        try %try load files number "nsave"
            velx=openvar('velx','NSTOKES_TEMP_vx_',nsave); % u velocity
            vely=openvar('vely','NSTOKES_TEMP_vy_',nsave); % v velocity
        catch %if error (files number "nsave" does not exist)
            %if calculations haulted in last section last saved set of frames
            %will sometimes be files number (nsave-1)
            velx=openvar('velx','NSTOKES_TEMP_vx_',nsave-1); % u velocity
            vely=openvar('vely','NSTOKES_TEMP_vy_',nsave-1);  % u velocity  
        end
        velx=velx(:,:,end); %get last timestep from external save u
        vely=vely(:,:,end); %get last timestep from external save v
        velx((bounds)==1)=nan; % set nodes inside solid wall region to nan (does not show up in velocity plot)
        vely((bounds)==1)=nan; % set nodes inside solid wall region to nan (does not show up in velocity plot)
        F=flipud(sqrt(velx.^2+vely.^2)); % velocity magnitude
        pcolor(x,y,F); %plot velocity magnitude for each node w.r.t x and y position
        xlabel(sprintf('x(meters)\n[time-step:%.0f of %.0f (last saved timestep)]',TLS,MI)) %label x axis
        ylabel('y(meters)') %label y axis
        title(sprintf('%s\nsimulation time=%.4fseconds',strmg1,TLS.*dt)) %title for figure (show time to 4 decimal places [%.4f])
        cb=colorbar; %create colorbar object with reference "cb"
        cba=caxis; %record current min/max values of colour bar (and hence min/max magnitude of velocity)
        ylabel(cb,'Velocity (m/s)') %add ylabel to colour bar
        if velfield==1 %if velocityfield is to be indicated (by streamslice() function)
            hd=streamslice(x,y,flipud(velx),-flipud(vely),arrowdensity); %plot streamslice (indicate velocity field)
            set(hd,'color',[0.0 0.0 0.0],'linewidth',arrowwidth) %set colour and line width of arrows in streamslice
        end
        sld=surf(x,y,ones(size(x)).*-1); %create foreground object (used to highlight solid boundaries with black colour)
        shading interp %reinforce smooth shading for pcolor and surf objects
        axis image %stretch axis to scale image
        set(sld,'facecolor',[0.0 0.0 0.0]) %set colour of sld to black
        caxis([cba(1) cba(2)]) %reset colorbar axis to only consider F (and not sld)
end
%###################################################################################################################################################



% PLOT OF VELOCITY AT SPECIFIC TIME STEP "TF"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if nsave~=0; %if external files were saved (nsave will not equal zero)
    %Determine time step to load and which external file to load#############################
        T2=1; %parameter which keeps track when next externally saved data (for u and v) should be loaded
        ns=1; %parameter to determine which file to load
        gotolast=0; %if =1, load very last frame;
        TF0=TF; %frame wanted to plot
        if TF>MI %if target frame exceeds what has been calculated
            gotolast=1; %load last frame with warning!
            TF=MI; %target frame is now last frame
        else %else if target frame does not exceed number of calculated time steps
            for tu=1:TF; %determine external file and frame to load
                T2=T2+1;   %update frame to load from current file 
                if T2>=ts %determine if next external file needs to be loaded
                        ns=ns+1;  %file number to load (next set of calculations)
                        T2=1+(T2-ts); %update time step parameter w.r.t to newly loaded files
                end
            end
            if ns>nsave; gotolast=1; end; %if ns greater than files saved in last run of section 1 (avoid using files from seperate calculations)
        end
        
    %######################################################################################
    figure(3) %plot of velocity for the last timestep saved in external file
        hold on %hold all objects created on figure (create new object without deleting previous)
        successload=0; %if =1, target frame was successfully loaded
        % LOAD NEEDED FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if gotolast==1; %if last frame is to be loaded [could not reach target frame]
            fprintf('For figure(3), frame %.0f was not calculated in section 1\n',TF0)
            fprintf('of the code (max number of frames is %.0f)\n',TLS)
            fprintf('so, final frame instead loaded (same as figure (2))\n')
            try %try load files number "nsave"
                velx=openvar('velx','NSTOKES_TEMP_vx_',nsave); % u velocity
                vely=openvar('vely','NSTOKES_TEMP_vy_',nsave); % v velocity
            catch %if error (files number "nsave" does not exist)
                %if calculations haulted in last section last saved set of frames
                %will sometimes be files number (nsave-1)
                velx=openvar('velx','NSTOKES_TEMP_vx_',nsave-1); % u velocity
                vely=openvar('vely','NSTOKES_TEMP_vy_',nsave-1);  % u velocity  
            end
            velx=velx(:,:,end); %get last timestep from external save u
            vely=vely(:,:,end); %get last timestep from external save v
            TLU=TLS; %frame that was loaded (target frame)
            successload=1; %frame successfully loaded (although it is same as figure (2))
        else
            try %try loaded needed external files and frame
                velx=openvar('velx','NSTOKES_TEMP_vx_',ns); % u velocity
                vely=openvar('vely','NSTOKES_TEMP_vy_',ns); % v velocity  
                velx=velx(:,:,T2); %get last timestep from external save u
                vely=vely(:,:,T2); %get last timestep from external save v
                successload=1; %file successfully loaded (without error above)
                TLU=TF; %frame that was loaded (last saved frame);
            catch %if could not load, flash warning in title of figure(3)
                str1=sprintf('Time-step %.0f could not be loaded as the externally saved files',TF);
                str2='could not be loaded or do not exist:';
                str3=sprintf('files; NSTOKES__TEMP__vx__%.0f.m and NSTOKES__TEMP__vy__%.0f.m',ns,ns);
                title(sprintf('%s\n%s\n%s',str1,str2,str3));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if successload==1; %if file was successfully loaded [plot figure]
            velx((bounds)==1)=nan; % set nodes inside solid wall region to nan (does not show up in velocity plot)
            vely((bounds)==1)=nan; % set nodes inside solid wall region to nan (does not show up in velocity plot)
            F=flipud(sqrt(velx.^2+vely.^2)); % velocity magnitude
            pcolor(x,y,F); %plot velocity magnitude for each node w.r.t x and y position
     
            xlabel(sprintf('x(meters)\n[time-step:%.0f of %.0f]',TLU,MI)) %label x axis
            ylabel('y(meters)') %label y axis
            title(sprintf('%s\nsimulation time=%.4fseconds',strmg1,TLU.*dt)) %title for figure (show time to 4 decimal places [%.4f])
            cb=colorbar; %create colorbar object with reference "cb"
            cba=caxis; %record current min/max values of colour bar (and hence min/max magnitude of velocity)
            ylabel(cb,'Velocity (m/s)') %add ylabel to colour bar
            if velfield==1 %if velocityfield is to be indicated (by streamslice() function)
                hd=streamslice(x,y,flipud(velx),-flipud(vely),arrowdensity); %plot streamslice (indicate velocity field)
                set(hd,'color',[0.0 0.0 0.0],'linewidth',arrowwidth) %set colour and line width of arrows in streamslice
            end
            sld=surf(x,y,ones(size(x)).*-10); %create foreground object (used to highlight solid boundaries with black colour)
            axis image %stretch axis to scale image
            shading interp %reinforce smooth shading for pcolor and surf objects
            set(sld,'facecolor',[0.0 0.0 0.0]) %set colour of sld to black
            caxis([cba(1) cba(2)]) %reset colorbar axis to only consider F (and not sld)
        end
end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@