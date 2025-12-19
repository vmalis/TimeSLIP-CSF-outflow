function listing = dir3(varargin)
%--------------------------------------------------------------------------
%
%   Subroutine to get list of folders with data
%
%   INPUT:  path to directory                                      'string'
%   OUTPUT: directory structure                                    struct
%__________________________________________________________________________
% VM (vmalis@ucsd.edu)
%--------------------------------------------------------------------------

if nargin == 0
    name = '.';
elseif nargin == 1
    name = varargin{1};
else
    error('Too many input arguments.');
end

% Get the directory listing
listing = dir(name);

% Filter out '.', '..', and '.DS_Store'
validEntries = ~ismember({listing.name}, {'.', '..', '.DS_Store'});

% Remove entries with invalid names
listing = listing(validEntries);

% Remove empty files but keep directories
nonEmptyOrDir = [listing.bytes] > 0 | [listing.isdir];
listing = listing(nonEmptyOrDir);

% Ensure only directories are listed
listing = listing([listing.isdir]);

end