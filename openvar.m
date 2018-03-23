function  G=openvar(var,str,n);
%Created by Jamie Johns 2016
% This function was created for use with the navier-stokes 2D incompressible flow code
% which can be found at;
%  https://github.com/JamieMJohns/Navier-stokes-2D-numerical-solve-incompressible-flow-with-custom-scenarios-MATLAB-

%The function of this file is to load previous calculated and saved
%x-component and y-component velocities (U and V).

%Inputs:
%   var=variable to load from temporary file
%   str=string which helps further specify which file to load
%   n=integer which helps further specify file to load
%Outputs:
%   G=variable which was stored in the loaded file

%note for below: pwd = current directory that matlab is "looking" at.
%another .m file (savevar.m, stores these files in directory "pwd"\temporary_NS_velocity)

    if exist([pwd '\temporary_NS_velocity'],'dir')==0; %if directory does not exits (for which saved files should be stored)
        G=nan; %no variable will be able to be loaded (set output to nan)
        fprintf('\nDirectory does not exist (related to openvar.m)\n') %prompt telling user that
    else %else, if the folder does exist
        g=load([pwd '\temporary_NS_velocity\' str num2str(n) '.mat' ],'var');%LOAD STORED VARIABLE (output is structure)
        G=g.var;%convert structure to variable (output of function file)   
    end    
    
end




