function mergedFolder = merge_dicom_TSLIP(rootDir)
if nargin==0, rootDir = pwd; end
D = dir(rootDir); D = D([D.isdir]); D = D(~ismember({D.name},{'.','..'}));
keep = ~cellfun('isempty', regexp({D.name}, '_(\d+)times')) & ~cellfun('isempty', regexp({D.name}, '_\d+$'));
D = D(keep); assert(~isempty(D),'No matching source folders.');

% Parse times & trailing suffix; sort by suffix
nt = cellfun(@(s) str2double(regexp(s,'_(\d+)times','tokens','once')), {D.name});
sx = cellfun(@(s) str2double(regexp(s,'_(\d+)$','tokens','once')), {D.name});
[~,ix] = sort(sx); D = D(ix); nt = nt(ix); sx = sx(ix);
sumT = sum(nt); firstS = sx(1); lastS = sx(end);

% Build MERGED folder name from FIRST folder base
nm = D(1).name;
bp = regexp(nm,'^[^_]+_(.*)_\d+times','tokens','once'); bp = bp{1};         % e.g., BRT_down_..._dT=500
Sp = regexp(nm,'(S\d+_10gap)','tokens','once'); if isempty(Sp), Sp = {''}; end; Sp = Sp{1};
mergedFolder = fullfile(rootDir, sprintf('MERGED_%d_%d_%s_%dtimes_%s', firstS, lastS, bp, sumT, Sp));
if ~exist(mergedFolder,'dir'), mkdir(mergedFolder); end

% Helper to list DICOM-like files (prefer *.dcm, else all non-dirs)
listImgs = @(p) local_list_imgs(p);

% 1) Copy FIRST folder as-is
L = listImgs(fullfile(rootDir, D(1).name)); assert(~isempty(L),'No images in first folder.');
for i = 1:numel(L)
    copyfile(fullfile(L(i).folder,L(i).name), fullfile(mergedFolder, L(i).name));
end

% 2) Derive naming pattern from the LAST copied file (to continue correctly)
lastName = L(end).name;
tok = regexp(lastName,'^(.*?)(\d+)(\.[^.]+)$','tokens','once');
if isempty(tok)  % fallback if no trailing number
    [~,bn,ex] = fileparts(lastName); tok = { [bn '_'], '0', iff(isempty(ex),'.dcm',ex) };
end
prefix = tok{1}; numWidth = numel(tok{2}); ext = tok{3}; counter = str2double(tok{2});

% 3) Copy from remaining folders with sequential renaming
for k = 2:numel(D)
    L = listImgs(fullfile(rootDir, D(k).name));
    for i = 1:numel(L)
        counter = counter + 1; % numeric increment (19 -> 20)
        newName = sprintf('%s%0*d%s', prefix, numWidth, counter, ext);
        copyfile(fullfile(L(i).folder,L(i).name), fullfile(mergedFolder,newName));
    end
end
fprintf('Merged %d folders -> %s\n', numel(D), mergedFolder);
end

% ------- local helpers -------
function L = local_list_imgs(p)
L = dir(fullfile(p,'*.dcm'));
if isempty(L), L = dir(p); L = L(~[L.isdir]); end
[~,ix] = sort({L.name}); L = L(ix);
end
function y = iff(c,a,b), if c, y=a; else, y=b; end
end