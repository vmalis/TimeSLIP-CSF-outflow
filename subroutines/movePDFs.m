function movePDFs(srcDir, dstDir)
%--------------------------------------------------------------------------
%
%   Function to move all PDF files from a source directory to a destination 
%   directory while preserving the directory structure.
%
%   INPUT:  
%       srcDir - 'string' 
%           Path to the source directory containing the PDF files.
%       dstDir - 'string'
%           Path to the destination directory where the PDF files will be 
%           copied to, with the original directory structure preserved.
%
%   OUTPUT: 
%       Moves all PDF files found in the source directory and its 
%       subdirectories to the destination directory, maintaining the 
%       folder structure. No explicit output is returned.
%
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%__________________________________________________________________________

    pdfFiles = dir(fullfile(srcDir, '**', '*.pdf'));

    % Loop over each PDF file found
    for k = 1:length(pdfFiles)
        % Get the relative path of the file from the source directory
        relativePath = strrep(pdfFiles(k).folder, srcDir, '');
        
        % Create the corresponding destination directory if it does not exist
        newDir = fullfile(dstDir, relativePath);
        if ~exist(newDir, 'dir')
            mkdir(newDir);
        end
        
        % Move the PDF file to the destination directory
        srcFile = fullfile(pdfFiles(k).folder, pdfFiles(k).name);
        dstFile = fullfile(newDir, pdfFiles(k).name);
        movefile(srcFile, dstFile);
    end

    disp('All PDFs were moved!');
end