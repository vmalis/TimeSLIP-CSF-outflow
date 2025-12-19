function [ARRAY,STRUCT] = dicom2struct_canon(path,filename)

%==========================================================================
% Subroutine to read dicom images to structure
%==========================================================================
%
% Read original images and create structure with Image, Location and header
% 
% INput:        path
% OUTput:       saves *.mat in path
%               outputs struct to workspace
%--------------------------------------------------------------------------
% written by Vadim Malis
% 04/16 at UCSD RIL
%==========================================================================
warning('off', 'images:dicominfo:fileVRDoesNotMatchDictionary')

cd(path)
dicomlist=dir('I*.dcm');
dicomlist=natsortfiles(dicomlist);
    
% get header info for current set from 1st image header
dicom_header = dicominfo(dicomlist(1).name);
numim=size(dicomlist,1);
r = dicom_header.Rows;
c = dicom_header.Columns;

h = waitbar(0,['loading image 1 out of ',num2str(numim)]);

for i=1:numim
	STRUCT(i).Image=dicomread(dicomlist(i).name);
	STRUCT(i).header=dicominfo(dicomlist(i).name);
%    STRUCT(i).ppm=STRUCT(i).header.ImageComments;
%     STRUCT(i).echoN = STRUCT(i).header.EchoNumber;
%     STRUCT(i).TE = STRUCT(i).header.EchoTime;
    STRUCT(i).location=STRUCT(i).header.SliceLocation;
    STRUCT(i).header = dicomCanonScaling(STRUCT(i).header);

    if isfield(STRUCT(i).header,'EchoNumber')
        STRUCT(i).echoN = STRUCT(i).header.EchoNumber;
        echos=max([STRUCT.echoN]);
    else
        echos=1;
    end
    
    if isfield(STRUCT(i).header,'EchoTime')
        STRUCT(i).TE = STRUCT(i).header.EchoTime;
    end
    
    if isfield(STRUCT(i).header,'SliceLocation')
        STRUCT(i).location=STRUCT(i).header.SliceLocation;
    end
    waitbar(i/numim,h,['loading image ', num2str(i), ' out of ', num2str(numim)])
end


echos=size(unique([STRUCT.TE]),2);
close(h)

ARRAY=reshape([STRUCT.Image],[r,c,numim/echos,echos]);
%ARRAY=reshape([STRUCT.Image],[r,c,numim]);

save(filename,'STRUCT')