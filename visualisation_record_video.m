%created by Jamie Johns 2018
%NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%This .m file is accessed and executed within "main.m" ( specifically "section4" after running
%"section 1").
%And running this .m file, independently, may not work; especially without first running the previously mentioned "sections" of main.m
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    WALLS=-5.*ones(size(x)); %matrix that will represtent solid walls
    WALLS(bounds==0)=nan; %if matrix element is not solid wall node (set to nan so is invisible)
    SOLIDWALL=surf(x,y,WALLS); %create surf object that will display solid wall nodes (coloured black)
    set(SOLIDWALL,'facecolor',[0,0,0]) %set surf object to black
    xlabel('X(meters)','fontsize',fntsize_axis_label) %set label for x axis
    ylabel('Y(meters)','fontsize',fntsize_axis_label) %set label for y axis
    view(0,-90) %flip plot view (vertically)
    xlim([0 domainX]) %set x-axis plot limits
    ylim([0 domainY]) %set y-axis plot limits
     if colourbar==1 %if colour bar to be displayed
        cba=colorbar; %create colour bar
        ylabel(cba,'Velocity magnitude (m/s)','fontsize',fntsize_axis_label_colorbar) %add y label to colour bar
      if capvel==1; caxis([0 mvel]); end; %if colour bar is kept constant w.r.t maximum calculated velocity (over all time steps)
     end
     axis image %set image axis to stretch to image specifications
title(sprintf('STRETCH FIGURE TO DESIRED RESOLUTION FOR RECORDED VIDEO\nThen press enter in matlab command window'))
fprintf('Press enter to begin recording animation\n')
fprintf('(or cntrl+c to cancle animation recording)\n')
user_prompt=input('<PRESS ENTER TO BEGIN RECORDING>','s');
%################################################################################################
     
     
%determine list of time-steps to plot u,v@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
num_frames_needed=timetarg.*framespersec; %number expected frames (running time x frames per second
    if MI>num_frames_needed; %if more timesteps were calculated than needed
        incf=round(MI./num_frames_needed); %determine rate to skip time-steps
        daframes=[1:incf:(MI-incf) MI];  %List of frames (time-steps) to plot
        repc=1; %number of times to repetitively record a frame (1 here since we are skipping frames)
    else %else, if number of time-steps calculated is less than needed
        daframes=1:MI; %list of all frames (timesteps)
        incf=1; %determine rate to skip frames (no skipping, go to next frame)
        repc=round(num_frames_needed/MI); %determine number of multipcle recordings for each frame
    end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


% CREATE OBJECT FOR WRITING VIDEO FILE$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if RECORD==1 %if animation is being recorded as video file
    vidnotexist=0; %parameter that determine if video of name does not exist
    vname=videoname; %initial video name
    parv=1; %suffix for video name, if previous file exists (i.e - vname_parv.avi)
    while vidnotexist==0;
        if exist([pwd '\' vname '.avi'])~=0 %check if video file of name vname already exist
            vname=[videoname '_' num2str(parv)]; %new file name with "_parv" (parv=integer)
            parv=parv+1; %update file number parameter
        else %else if file does not exist
            vidnotexist=1; % exit while loop (video filename is unique in folder)
        end
    end
    video_object = VideoWriter([vname '.avi']); %create object to write recorded frames to
    video_object.FrameRate = framespersec; %set frame rate of video file
    video_object.Quality = 100;   %set quality of video (100=100 percent quality)
    open(video_object);  %open video object to write frames to
    %Dialogue describing video file%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\n A video (avi) file will be created from the animation:\n')
        fprintf('->Target run time length: %.2f seconds (actual time: %.2fseconds)\n',timetarg,repc.*length(daframes)./framespersec)
        fprintf('->frame rate:%.2f frame per second\n',framespersec)
        fprintf('\nVideo file is named:\n%s\n',[vname '.avi'])
        fprintf('\nAnd will saved in folder:\n%s\n\n',pwd)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else %else if video is not being recorded from animated figure
    fprintf('\nNo video is being recorded from animated figure\n')
    fprintf('(parameter "Record" was set to 0 by the user)\n')
end
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


T2=1; %parameter which keeps track when next externally saved data (for u and v) should be loaded
ns=1; %parameter which keeps track which externally saved data to load
%Below is part of string used in title in figure (STATIC PART OF STRING)
strmg1=sprintf('Rendered in Matlab with code created by Jamie M. Johns\n density=%.2fkg/m^3; \\mu=%.4fkg/(m*s); dt=%.4fs, resolution:%.0fx%.0f [for calculations]',dens,mu,dt,XI,YI);
title(sprintf('%s',strmg1),'fontsize',fntsize_title) %initial title for animation

vel_field=[]; %initialise parameter that will be object for plotted streamslice (velocity field)
M=sqrt(velx.^2+vely.^2); %velocity magnitude

for T=1:(length(daframes)); %for each time step to be plotted

    set(vel_mag,'cdata',M(:,:,T2)); %update colour data to current velocity magnitudes

    if velfield==1; %if indicate velocityfield (using stream slice)$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if ~isempty(vel_field); delete(vel_field); clear hd; end; %delete vel_field (from last time-step)
        vel_field=streamslice(x,y,velx(:,:,T2),vely(:,:,T2),arrowdensity); %streamslice and set arrow/line density
       set(vel_field,'color',[0.0 0.0 0.0],'linewidth',arrowwidth) %set colour (to black) and width of lines
    end%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    %update title to current time step information (strmg1=static string in title):
    if title_detailed==1; %if using detailed title
        title(sprintf('%s\n time=%.4f seconds; frame:%.0f/%.0f',strmg1,daframes(T).*dt,T,length(daframes)),'fontsize',fntsize_title)
    end
    drawnow %update plot objects visually on figure
    
    if RECORD==1 %if video being recorded
        current_frame = getframe(figure(1)); %capture frame
        for rep=1:repc; %record same frame "repc" number of times
            writeVideo(video_object,current_frame); %write frame to video record object
        end;
    end;
    
    T2=T2+incf; %update parameter which keeps track of next time step w.r.t to current loaded file
    if T2>=ts %determine if next external file needs to be loaded
        ns=ns+1;  %file number to load for next set of calculations
        if ns<=nsave %if not at maximum number of saved files
           velx=openvar('velx','NSTOKES_TEMP_vx_',ns); %load next set of saved u
           vely=openvar('vely','NSTOKES_TEMP_vy_',ns); %load next set of saved v
           M=sqrt(velx.^2+vely.^2); %calculate velocity magnitude (for next set of u an v)
        end
       T2=1+(T2-ts); %update time step parameter w.r.t to newly loaded files
       if T==length(daframes); T2=ts; end; %if last frame make sure T refers very last frame of very last file
    end

end 

if RECORD==1 %if video is being recorded
    close(video_object); %close writer object
    fprintf('\n Video Recording has finished!\n')
    else %else, if video not being recorded
    fprintf('\n Animation has finished!\n')
end  