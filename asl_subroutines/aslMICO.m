function ROI = aslMICO(seriesPath, param, doplot)
% ROI = roi_from_T1w(catPath, param, doplot)
% - Looks for catPath/T1w/<series with DICOMs>
% - Runs brainMICO and returns an ROI struct with:
%     ROI.maskVm : [X,Y,3]  (2D) or [X,Y,Z,3] (3D)  -> WM,GM,CSF along last dim
%     ROI.maskIm : [X,Y,3]  single slice (WM,GM,CSF) at chosen slice_number
%     ROI.slice_number : chosen slice index (see logic below)
%     ROI.name   : {'WM','GM','CSF'}
%     ROI.area   : NaN  (placeholder)
%
% Return [] if no T1w or no DICOMs found.


ROI = [];

volPath = seriesPath;

% --- read volume
[V,~] = read_dicom_volume(volPath);          % V: XxYxZ (Z=1 for 2D)
[nr,nc,nz] = size(V);
V=int16(V);

maskall=logical(MaskInteractive(V));
V=V.*int16(maskall);

% default param.maskall if not provided/mismatched
if isempty(param) || ~isfield(param,'maskall') || ~isequal(size(param.maskall), size(V))
    param.maskall = true(size(V));
end

% --- run your segmentation
[mask_wm, mask_gm, mask_csf] = brainMICO(V, param, doplot);  % XxYxZ logical

% --- build maskVm (stack labels on last dim)
if nz == 1
    maskVm = cat(3, double(mask_wm), double(mask_gm), double(mask_csf));      % [X,Y,3]
else
    maskVm = cat(4, double(mask_wm), double(mask_gm), double(mask_csf));      % [X,Y,Z,3]
end

% --- choose slice_number for maskIm:
% first slice that has any mask; else middle slice
anyMask = squeeze(any(any(mask_wm|mask_gm|mask_csf,1),2));  % 1xZ
idx = find(anyMask,1,'first'); if isempty(idx), idx = ceil(nz/2); end
sl = idx;

% --- build maskIm (single slice, labels along 3rd dim)
maskIm = cat(3, double(mask_wm(:,:,sl)), double(mask_gm(:,:,sl)), double(mask_csf(:,:,sl))); % [X,Y,3]

% --- assemble ROI to return
ROI = struct();
ROI.maskVm = maskVm;            % [X,Y,3] or [X,Y,Z,3]
ROI.maskIm = maskIm;            % [X,Y,3]
ROI.slice_number = sl;
ROI.name = {'WM','GM','CSF'};
ROI.area = NaN;                 % placeholder to match your struct

end

% ---------- helper ----------
function [V,info0] = read_dicom_volume(folder)
L = dir(fullfile(folder,'*.dcm'));
assert(~isempty(L),'No DICOMs in %s', folder);
[~,ix] = sort({L.name}); L = L(ix);
p0 = fullfile(L(1).folder,L(1).name);
I0 = dicomread(p0); info0 = dicominfo(p0); %#ok<NASGU>
nr = size(I0,1); nc = size(I0,2); nz = numel(L);
V = zeros(nr,nc,nz, class(I0));
V(:,:,1) = I0;
for k = 2:nz
    V(:,:,k) = dicomread(fullfile(L(k).folder,L(k).name));
end
end