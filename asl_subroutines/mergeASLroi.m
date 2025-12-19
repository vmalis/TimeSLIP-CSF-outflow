function ROI = mergeASLroi(ROI, R)
% MERGEASLROI  Append ROI R (WM/GM/CSF) into ROI, resizing XY only.
% - ROI.maskVm: [X,Y,L]  or [X,Y,Z,L]   (labels last)
% - ROI.maskIm: [X,Y,L]
% - R.maskVm  : [x,y,3]  or [x,y,z,3]
% - R.maskIm  : [x,y,3]
% - R.name    : {'WM','GM','CSF'}

if isempty(R) || ~isfield(R,'maskVm') || isempty(R.maskVm), return; end

% If no existing ROI, just adopt R (keep storage labels-last)
if ~exist('ROI','var') || isempty(ROI) || ~isfield(ROI,'maskVm') || isempty(ROI.maskVm)
    ROI = R;
    return
end

% ---- Target geometry from existing ROI ----
tgtX = size(ROI.maskVm,1);
tgtY = size(ROI.maskVm,2);
tgtIs4D = (ndims(ROI.maskVm) >= 4);             % true if ROI is [X,Y,Z,L]
if tgtIs4D
    tgtZ = size(ROI.maskVm,3);
else
    tgtZ = 1;                                   % treat 2D as Z=1
end

% ---- Prepare R.maskVm to match ROI's XY (and Z length if ROI is 4D) ----
Rvm = R.maskVm;
if ndims(Rvm)==3                 % [x,y,3]  (2D labels)
    % resize XY per label to [tgtX,tgtY]
    out = zeros(tgtX, tgtY, size(Rvm,3), 'like', ROI.maskVm);
    for l = 1:size(Rvm,3)
        out(:,:,l) = imresize(Rvm(:,:,l), [tgtX tgtY], 'nearest') > 0.5;
    end
    if tgtIs4D
        % expand to Z=1 for concatenation with [X,Y,Z,L]
        Rvm_ready = reshape(out, tgtX, tgtY, 1, size(out,3));
    else
        Rvm_ready = out;  % [X,Y,3]
    end

elseif ndims(Rvm)>=4             % [x,y,z,3]
    rz = size(Rvm,3);
    % resize XY per (z,label)
    tmp = zeros(tgtX, tgtY, rz, size(Rvm,4), 'like', ROI.maskVm);
    for l = 1:size(Rvm,4)
        for z = 1:rz
            tmp(:,:,z,l) = imresize(Rvm(:,:,z,l), [tgtX tgtY], 'nearest') > 0.5;
        end
    end
    if tgtIs4D
        % match Z length by clip/pad (repeat last)
        if rz > tgtZ
            Rvm_ready = tmp(:,:,1:tgtZ,:);
        elseif rz < tgtZ
            pad = repmat(tmp(:,:,end,:), 1,1, tgtZ-rz, 1);
            Rvm_ready = cat(3, tmp, pad);
        else
            Rvm_ready = tmp;  % same Z
        end
    else
        % ROI is 2D -> pick one slice
        if isfield(ROI,'slice_number') && ~isempty(ROI.slice_number)
            zpick = min(max(1, round(ROI.slice_number)), rz);
        else
            % choose first slice that has any mask; fallback to 1
            anyMaskZ = squeeze(any(any(tmp,1),2));   % [z,l]
            zpick = find(any(anyMaskZ,2),1,'first'); if isempty(zpick), zpick = 1; end
        end
        Rvm_ready = squeeze(tmp(:,:,zpick,:));       % [X,Y,3]
    end

else                              % unexpected [x,y] -> promote to 3 labels
    base = imresize(Rvm, [tgtX tgtY], 'nearest') > 0.5;
    out = cat(3, base, base, base);
    if tgtIs4D, Rvm_ready = reshape(out, tgtX, tgtY, 1, 3);
    else,       Rvm_ready = out; end
end

% ---- Resize R.maskIm to ROI.maskIm XY and stack along 3rd dim ----
if ~isfield(ROI,'maskIm') || isempty(ROI.maskIm)
    imX = tgtX; imY = tgtY;
else
    imX = size(ROI.maskIm,1); imY = size(ROI.maskIm,2);
end
Rmi = R.maskIm;
if isempty(Rmi)
    Rmi_ready = zeros(imX, imY, 3, 'like', ROI.maskIm);
else
    Rmi_ready = zeros(imX, imY, size(Rmi,3), 'like', ROI.maskIm);
    for l = 1:size(Rmi,3)
        Rmi_ready(:,:,l) = imresize(Rmi(:,:,l), [imX imY], 'nearest') > 0.5;
    end
end

% ---- Concatenate along label dim ----
if tgtIs4D
    % ROI: [X,Y,Z,L] ; Rvm_ready: [X,Y,Z,3]
    ROI.maskVm = ROI.maskVm + sum(Rvm_ready,4);
else
    % ROI: [X,Y,L] ; Rvm_ready: [X,Y,3]
    ROI.maskVm = ROI.maskVm + sum(Rvm_ready,3);
end
    ROI.maskVm(ROI.maskVm>1)=1;

if ~isfield(ROI,'maskIm') || isempty(ROI.maskIm)
    ROI.maskIm = Rmi_ready;                % initialize
else
    ROI.maskIm = cat(3, ROI.maskIm, Rmi_ready);
end

% ---- Names ----
if ~isfield(ROI,'name') || isempty(ROI.name), ROI.name = {}; end
ROI.name = [ROI.name(:).'  R.name(:).'];

end