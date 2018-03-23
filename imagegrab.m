function [RED BLACK GREEN BLUE XI YI]=imagegrab(imageinput,xinc);
%Created by Jamie Johns 2016
% This file was made for use with Navier-stokes 2d incompressible flow
% code, which can be found at: (github link)
% https://github.com/JamieMJohns/Navier-stokes-2D-numerical-solve-incompressible-flow-with-custom-scenarios-MATLAB-

%Note: matlab must be "looking" at same directory as used image for "imageinput"
%      as well as for this function file (put this function file and image to be used in same directory).

%The purpose of this code is to extract details from input images;
%Inputs:
%   imagegrap = input image (e.g - imageinput='Myphoto.jpg")
%   xinc= width resolution of mask created from image ( note- output XI is same as YI)
%   -> output YI is height resolution of output masks which is scale
%      with respect to dimensions of input image and specified "xinc"
%    e.g - if, imagegrab='myphoto.png' has resolution of 1000x500 (WxL; W/L=2)
%                  xinc=200
%          then, XI=200 and YI=100
%          therefore, size(RED)=size(BLACK)=size(BLUE)=200x100  (masks created for downsampled image [downsample image created from original image])
%Outputs;
%   RED=binary map mask for input image, RED(i,j)=1 if i,j corresponds to red colored pixel ([R G B]~=[255 0 0]) on downsampled image (and =0 if not)
%   BLACK=binary map mask for input image, BLACK(i,j)=1 if i,j corresponds to black colored pixel ([R G B]~=[0 0 0]) on downsampled image (and =0 if not)
%   GREEN=binary map mask for input image, GREEN(i,j)=1 if i,j corresponds to green colored pixel ([R G B]~=[0 255 0]) on downsampled image (and =0 if not)
%   BLUE=binary map mask for input image, BLUE(i,j)=1 if i,j corresponds to blue colored pixel ([R G B]~=[0 255 0]) on downsampled image (and =0 if not)
%   XI=width resolution for masks (and hence width resolution of downsampled [or "upsampled"] image, imageinput)
%   YI=height resolution for masks (and hence width resolution of downsampled [or "upsampled"] image, imageinput)



%RENAME VARIABLE HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
im= imread(imageinput); %load colormap of image ' '(must be looking in!!!!!!! im=3d matrix; pixel intensity (R,G,B) for each pixel of input image
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
imdim=size(im); %original resolution of input image
yinc=round((imdim(1)./imdim(2)).*xinc); %determine new vertical resolution of image
xl=1:imdim(2); %horizontal indexing for input image
yl=1:imdim(1); %vertical indexing for input image
xj=linspace(1,imdim(2),xinc); %horizontal indexing for downsampled image
yj=linspace(1,imdim(1),yinc); %vertical indexing for downsampled image
XI=xinc; %number of horizontal pixels in downsampled image
YI=yinc; %number of vertical pixels in downsampled image
imnw=zeros(yinc,xinc,3); %pixels for downsampled imaged

% For downsampled image pixel, use pixel from original image that is closest i proximity
for j=1:xinc; %for each downsampled horizontal pixel
    for k=1:yinc; %for each downsampled vertical pixe
     [a bx]=min(abs(xl-xj(j))); %find closest horizontal pixel in original image
     [a by]=min(abs(yl-yj(k))); %find closest vertical pixel in original image
     imnw(k,j,1)=im(by,bx,1); %record closest proximity (by "distance") for color red
      imnw(k,j,2)=im(by,bx,2); %record closest proximity (by "distance") for color green
      imnw(k,j,3)=im(by,bx,3); %record closest proximity (by "distance") for color blue
    end
end



on=ones(yinc,xinc); %matrix of ones with same resolution as downsampled (DS) imaged [kept in use below as placeholder for further operations(if edit happens)]
k=1:XI; %index of horizontal pixels in downsampled image
j=1:YI; %index of vertical pixels in downsampled image
im=double(im); %convert original matrix to double precision
sharpen=1; % if use "sharpening" (below), to sharpen white and black in image

% saturation of colors############################# [note: 0=nothing (mininum intensity) 255=maximum intensity)
%apply thresholding (saturation) to further
            %Saturate yellows (if R,G>200 and B<150 ; set R,G=255 and B=0) [note yellow masks are not currently created from this function file]
            testr=(imnw(j,k,1)>=200).*(imnw(j,k,2)>=200).*(imnw(j,k,3)<=150); % =1 , if apply saturation
             imnw(j,k,1)=testr.*255.*on+(~testr).*imnw(j,k,1); %operations to apply to red
             imnw(j,k,2)=testr.*255.*on+(~testr).*imnw(j,k,2); %operations to apply to green
             imnw(j,k,3)=(~testr).*imnw(j,k,3); %operations to apply to blue
             %Saturate Green (if G>200 and R,B<150 ; set G=255 and R,B=0) 
            testr(j,k)=(imnw(j,k,1)<=150).*(imnw(j,k,2)>=200).*(imnw(j,k,3)<=150); % =1 , if apply saturation
            imnw(j,k,1)=(~testr(j,k)).*imnw(j,k,1); %operations to apply to red
            imnw(j,k,2)=testr.*255+(~testr).*imnw(j,k,2); %operations to apply to green
            imnw(j,k,3)=(~testr).*imnw(j,k,3); %operations to apply to blue
            %Saturate Blue (if B>200 and R,G<150 ; set B=255 and R,G=0) 
            testr(j,k)=(imnw(j,k,1)<=150).*(imnw(j,k,3)>=200).*(imnw(j,k,2)<=150); % =1 , if apply saturation
            imnw(j,k,1)=(~testr(j,k)).*imnw(j,k,1); %operations to apply to red
            imnw(j,k,3)=testr.*255+(~testr).*imnw(j,k,3); %operations to apply to blue
            imnw(j,k,2)=(~testr).*imnw(j,k,2);%operations to apply to green
            %Saturate Red (if R>200 and B,G<150 ; set R=255 and B,G=0) 
            testr(j,k)=(imnw(j,k,1)>=200).*(imnw(j,k,2)<=150).*(imnw(j,k,3)<=150); % =1 , if apply saturation
            imnw(j,k,2)=(~testr(j,k)).*imnw(j,k,2);%operations to apply to green
            imnw(j,k,1)=testr.*255+(~testr).*imnw(j,k,1); %operations to apply to red
            imnw(j,k,3)=(~testr).*imnw(j,k,3); %operations to apply to blue
            %Saturate cyans (if B,G>200 and R<150 ; set B,G=255 and R=0) [note cyan masks are not currently created from this function file]
             testr=(imnw(:,:,1)<=150).*(imnw(:,:,2)>=200).*(imnw(:,:,3)>=200); % =1 , if apply saturation
             imnw(:,:,3)=testr.*255.*on+(~testr).*imnw(:,:,3); %operations to apply to blue
             imnw(:,:,2)=testr.*255.*on+(~testr).*imnw(:,:,2); %operations to apply to green
             imnw(:,:,1)=(~testr).*imnw(:,:,1); %operations to apply to red
% sharpen colours (define of white and black)#########################
             if sharpen~=0 %if sharpening between white and black to be applied
              % if (R,G,B)>150 ,set (R,G,B)=255 [completely white in colour]
              testr=(imnw(:,:,1)>=150).*(imnw(:,:,2)>=150).*(imnw(:,:,3)>=150); % =1 , if apply "sharpening"
             imnw(:,:,3)=testr.*250.*on+(~testr).*imnw(:,:,3);%operations to apply to blue
             imnw(:,:,2)=testr.*250.*on+(~testr).*imnw(:,:,2);%operations to apply to green
             imnw(:,:,1)=testr.*250.*on+(~testr).*imnw(:,:,1);
             % if (R,G,B)<150 ,set (R,G,B)=0 [completely black in colour]
             testr=(imnw(:,:,1)<=150).*(imnw(:,:,2)<=150).*(imnw(:,:,3)<=150);% =1 , if apply "sharpening"
             imnw(:,:,3)=(~testr).*imnw(:,:,3);%operations to apply to blue
             imnw(:,:,2)=(~testr).*imnw(:,:,2);%operations to apply to green
             imnw(:,:,1)=(~testr).*imnw(:,:,1);%operations to apply to red
             end


%convert imnw data to 2d per-colour masks (ouput of this function file)
  RED=(imnw(:,:,2)==0).*(imnw(:,:,3)==0).*(imnw(:,:,1)>=200); %Red mask for downsampled image (RED(i,j)=1 if pixel i,j is red colour)
  GREEN=(imnw(:,:,1)==0).*(imnw(:,:,2)>=150).*(imnw(:,:,3)==0);%Green mask for downsampled image (RED(i,j)=1 if pixel i,j is red colour)
  BLACK=(imnw(:,:,2)==0).*(imnw(:,:,3)==0).*(imnw(:,:,1)==0);%Black mask for downsampled image (BLACK(i,j)=1 if pixel i,j is black colour)
  BLUE=(imnw(:,:,2)==0).*(imnw(:,:,3)>=150).*(imnw(:,:,1)==0);%Blue mask for downsampled image (BLUE(i,j)=1 if pixel i,j is blue colour)
           
            
end           
            