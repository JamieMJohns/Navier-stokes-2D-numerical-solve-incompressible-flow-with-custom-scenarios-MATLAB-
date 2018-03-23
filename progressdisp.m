function vt=progressdisp(i,maxi,str,str2,vi,time);
%Created by Jamie Johns 2016
%This code is created for use of my Navier stokes 2d incompressible flow
%which can be found at;

%https://github.com/JamieMJohns/Navier-stokes-2D-numerical-solve-incompressible-flow-with-custom-scenarios-MATLAB-

%In relation to the files which can be found in the link above;
%this function file serves as easy to use (and re-use)
%"progress of calculations" dialogue in matlab command window
% (increments of 10%) and to also created to save space (lines) in
% main.m.


%Inputs;
%   i=paramter that keeps track of current instant (ith instant)
%   maxi=maximum number of instances
%   str=string for dialogue (example given below)
%   str2=another string for dialogue (example given below)
%   vi=parameter which keeps track of which string to use below
%   time=time in seconds (typical "toc" is used in time);
%Outputs;
%   vt=update of dialogue parameter;
%   string printed in command window, if condition satisfied.

%e.g- vt=progressdisp(44,200,'calculations','processed',2,13.42);
%     ouput;
%       Printed in command window;
%           "20.00% of calculations completed at 13.4200 seconds [44 of 200 total processed]"   
%       [44 is 22% of  200 (dialogue [above] just indicates that at least 20% complete)]
%       vt=vi+1=2+1=3

%NOTE: the first time you run this file (such as in for loop); set vt=0


    vt=vi; %initialize vt (value of vt if condition not satisfied below
    if vt==0 %if first time running progressdisp() [typically vt=vi=0]
        fprintf('\n\n %s started at %.4fseconds',str,time) %initial file use dialogue
        vt=1; %update dialogue parameter (when vt is now used as input "vi" [when function next used], 10% complete dialogue is next to be printed)
    end    

    if i>=((vt./10)*maxi) && vt==1; %if at least 10% of instances processed; (i/maxi)>=0.1 and vt=1
        fprintf('\n %.2f%% of %s completed at %.4fseconds [%.0f of %.0f total %s]',(vt./10).*100,str,time,i,maxi,str2) %print 10% COMPLETE dialogue
        vt=2;    %update dialogue parameter (when vt is now used as input "vi" [when function next used], 20% complete dialogue is next to be printed)
    elseif i>=((vt./10)*maxi) && vt==2; %if at least 20% of instances processed; (i/maxi)>=0.2 and vt=2
        fprintf('\n %.2f%% of %s completed at %.4fseconds [%.0f of %.0f total %s]',(vt./10).*100,str,time,i,maxi,str2) %print 20% COMPLETE dialogue
        vt=3;    %update dialogue parameter (when vt is now used as input "vi" [when function next used], 30% complete dialogue is next to be printed)   
    elseif i>=((vt./10)*maxi) && vt==3; %if at least 30% of instances processed; (i/maxi)>=0.3 and vt=3
        fprintf('\n %.2f%% of %s completed at %.4fseconds [%.0f of %.0f total %s]',(vt./10).*100,str,time,i,maxi,str2) %print 30% COMPLETE dialogue
        vt=4;    %update dialogue parameter (when vt is now used as input "vi" [when function next used], 40% complete dialogue is next to be printed)
    elseif i>=((vt./10)*maxi) && vt==4; %if at least 40% of instances processed; (i/maxi)>=0.4 and vt=4
        fprintf('\n %.2f%% of %s completed at %.4fseconds [%.0f of %.0f total %s]',(vt./10).*100,str,time,i,maxi,str2) %print 40% COMPLETE dialogue
        vt=5;    %update dialogue parameter (when vt is now used as input "vi" [when function next used], 50% complete dialogue is next to be printed)
    elseif i>=((vt./10)*maxi) && vt==5; %if at least 50% of instances processed; (i/maxi)>=0.5 and vt=5
        fprintf('\n %.2f%% of %s completed at %.4fseconds [%.0f of %.0f total %s]',(vt./10).*100,str,time,i,maxi,str2) %print 50% COMPLETE dialogue
        vt=6;    %update dialogue parameter (when vt is now used as input "vi" [when function next used], 60% complete dialogue is next to be printed)
    elseif i>=((vt./10)*maxi) && vt==6; %if at least 60% of instances processed; (i/maxi)>=0.6 and vt=6
        fprintf('\n %.2f%% of %s completed at %.4fseconds [%.0f of %.0f total %s]',(vt./10).*100,str,time,i,maxi,str2) %print 60% COMPLETE dialogue
        vt=7;    %update dialogue parameter (when vt is now used as input "vi" [when function next used], 70% complete dialogue is next to be printed)
    elseif i>=((vt./10)*maxi) && vt==7; %if at least 70% of instances processed; (i/maxi)>=0.7 and vt=7
        fprintf('\n %.2f%% of %s completed at %.4fseconds [%.0f of %.0f total %s]',(vt./10).*100,str,time,i,maxi,str2) %print 70% COMPLETE dialogue
        vt=8;    %update dialogue parameter (when vt is now used as input "vi" [when function next used], 80% complete dialogue is next to be printed)
    elseif i>=((vt./10)*maxi) && vt==8; %if at least 80% of instances processed; (i/maxi)>=0.8 and vt=8
        fprintf('\n %.2f%% of %s completed at %.4fseconds [%.0f of %.0f total %s]',(vt./10).*100,str,time,i,maxi,str2) %print 80% COMPLETE dialogue
        vt=9;    %update dialogue parameter (when vt is now used as input "vi" [when function next used], 90% complete dialogue is next to be printed)
    elseif i>=((vt./10)*maxi) && vt==9; %if at least 90% of instances processed; (i/maxi)>=0.9 and vt=9
        fprintf('\n %.2f%% of %s completed at %.4fseconds [%.0f of %.0f total %s]',(vt./10).*100,str,time,i,maxi,str2) %print 90% COMPLETE dialogue
        vt=10;    %update dialogue parameter (when vt is now used as input "vi" [when function next used], 100% complete dialogue is next to be printed)
    elseif i>=((vt./10)*maxi) && vt==10; %if at least 100% of instances processed; (i/maxi)>=1 and vt=10
        fprintf('\n %.2f%% of %s completed at %.4fseconds [%.0f of %.0f total %s]',(vt./10).*100,str,time,i,maxi,str2) %print 100% COMPLETE dialogue
        vt=11;    %update dialogue parameter (when vt is now used as input "vi" [when function next used]; no further dialogue will be printed)
    else
        %no other options (place keeper)
    end



end