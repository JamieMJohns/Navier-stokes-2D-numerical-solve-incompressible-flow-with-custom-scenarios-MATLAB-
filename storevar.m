function  G=storevar(var,str,n)
%Created by Jamie Johns 2016
% This function was created for use with the navier-stokes 2D incompressible flow code
% which can be found at;

%https://github.com/JamieMJohns/Navier-stokes-2D-numerical-solve-incompressible-flow-with-custom-scenarios-MATLAB-

%The function of this file is to save calculated
%x-component and y-component velocities (U and V), for later usage.

%Inputs:
%   var=variable to save to file
%   str=string which helps further specify name of file that is saved
%   n=integer which helps further specify name of file that is saved
%Outputs:
%   No outputs (G is unused)

%note for below: pwd = current directory that matlab is "looking" at.
if exist([pwd '\temporary_NS_velocity'],'dir')==0; %if directory does not exits (for which saved files should be stored)
    mkdir([pwd '\temporary_NS_velocity']) %create desired directory for which files will be saved to.
end  
save([pwd '\temporary_NS_velocity\' str num2str(n) '.mat' ],'var') %save variable to file of containing particular string ("str") and integer ("int")
            
end




