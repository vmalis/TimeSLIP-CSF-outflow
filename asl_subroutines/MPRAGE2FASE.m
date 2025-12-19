function matched_slice = MPRAGE2FASE(mprage_folder, fase2d_folder)
%==========================================================================
% MPRAGE2FASE
%--------------------------------------------------------------------------
% Extract a slice from MPRAGE that matches the single-slice FASE2D geometry.
% Uses: dicom2struct_canon, dicomVolumeGeometry, affine, affine2vert,
%       dicomPlotVolumeOrientation, export_fig (unchanged).
% Creates Vol_Orientation.pdf and Slice_Match.png for QA.
%==========================================================================

%% Read MPRAGE
cwd = pwd;
cd(mprage_folder);
[volM, infoM] = dicom2struct_canon(pwd,'data');   % 3D
cd(cwd);

%% Read FASE2D
cd(fase2d_folder);
[imgF, infoF] = dicom2struct_canon(pwd,'data');   % single-slice series
cd(cwd);

% Ensure FASE is 3D (nz=1) for consistent 3D math
if ndims(imgF)==2
    volF = reshape(imgF, size(imgF,1), size(imgF,2), 1);
else
    volF = imgF;
end

% Normalize (same style as your code)
volM = mat2gray(volM);
volF = mat2gray(volF);

%% Geometry (keep dicomVolumeGeometry unchanged; only ask needed outputs)
[~, AM] = dicomVolumeGeometry(infoM);   % 4x4
[~, AF] = dicomVolumeGeometry(infoF);   % 4x4

%% Spacings from headers (for MPRAGE resample target)
dxM = double(infoM(1).header.PixelSpacing(1));
dyM = double(infoM(1).header.PixelSpacing(2));
if isfield(infoM(1).header,'SpacingBetweenSlices')
    dzM = double(infoM(1).header.SpacingBetweenSlices);
elseif isfield(infoM(1).header,'SliceThickness')
    dzM = double(infoM(1).header.SliceThickness);
else
    dzM = [];
end
if isempty(dzM) && numel(infoM)>=2
    T1 = double(infoM(1).header.ImagePositionPatient(:));
    T2 = double(infoM(2).header.ImagePositionPatient(:));
    io = double(infoM(1).header.ImageOrientationPatient(:));
    ez = cross(io(1:3), io(4:6));
    dzM = abs(dot(T2-T1, ez));
end
if isempty(dzM), dzM = 1.0; end

% FASE in-plane spacings (no resample of FASE via affine)
dxF = double(infoF(1).header.PixelSpacing(1));
dyF = double(infoF(1).header.PixelSpacing(2));

% Target spacing for MPRAGE (min XY across series, keep MPRAGE Z)
delta = [min(dxM,dxF), min(dyM,dyF), dzM];

%% Resample MPRAGE only (3D)
[volM_r, AM_r] = affine(volM, AM, delta, 2);

%% -------- ORIENTATION PDF (MPRAGE box + FASE 4-vertex plane) -----------
% MPRAGE (resampled) vertices in world coords (8 corners)
[vXM, vYM, vZM] = affine2vert(AM_r, size(volM_r));

% FASE plane: compute the 4 in-plane corners by mapping 0-based indices
nR = size(volF,1); 
nC = size(volF,2);
corners_idx = [ ...
    0      nR-1   nR-1   0   ;  % rows (i)
    0      0      nC-1   nC-1;  % cols (j)
    0      0      0      0   ;  % single slice (k=0)
    1      1      1      1   ]; % homogeneous
corners_xyz = AF * corners_idx;                 % 4x4
X4 = corners_xyz(1,:); 
Y4 = corners_xyz(2,:); 
Z4 = corners_xyz(3,:);

% Build a padded "full volume" box for context
d = 1;
XoMIN = min([vXM, X4]) - d;  XoMAX = max([vXM, X4]) + d;
YoMIN = min([vYM, Y4]) - d;  YoMAX = max([vYM, Y4]) + d;
ZoMIN = min([vZM, Z4]) - d;  ZoMAX = max([vZM, Z4]) + d;

Xbox = repmat([XoMIN, XoMAX],1,4);
Ybox = repmat([YoMIN, YoMIN, YoMAX, YoMAX],1,2);
Zbox = [repmat(ZoMIN,1,4), repmat(ZoMAX,1,4)];

% ---- Plot (same aesthetic as volumeAlign), FASE as a single patch ----
colors = distinguishable_colors(3);

figure('Position',[0,0,800,600]); set(gcf,'color','w'); view(3);
axis vis3d; axis square;
xlim([XoMIN*1.2, XoMAX*1.2]); ylim([YoMIN*1.2, YoMAX*1.2]); zlim([ZoMIN*1.2, ZoMAX*1.2]);
hold on; grid on;

hAxes = gca; hAxes.TickLabelInterpreter = 'latex'; hAxes.FontSize = 12;
xlabel('sagittal [mm]','FontSize',16,'Interpreter','latex');
ylabel('coronal [mm]','FontSize',16,'Interpreter','latex');
zlabel('axial [mm]','FontSize',16,'Interpreter','latex');

% MPRAGE (resampled) as a box
ph(1) = dicomPlotVolumeOrientation(vXM, vYM, vZM, colors(1,:), 0.40);

% FASE as a true 4-vertex plane
ph(2) = patch('XData',X4,'YData',Y4,'ZData',Z4, ...
              'FaceColor',colors(2,:), 'FaceAlpha',0.50, ...
              'EdgeColor',colors(2,:), 'LineWidth',1.5);

% Full volume context box
ph(3) = dicomPlotVolumeOrientation(Xbox, Ybox, Zbox, colors(3,:), 0.00);

legend(ph, {'MPRAGE (resampled)','FASE2D plane','Full volume'}, ...
       'FontSize',14,'Interpreter','latex');

export_fig('Vol_Orientation.pdf','-dpdf','-nocrop','-opengl');
close;

%% Build FASE pixel grid (0-based, k=0 plane), map → world → MPRAGE index
nR = size(volF,1);  nC = size(volF,2);
[RR0, CC0] = ndgrid(0:nR-1, 0:nC-1);
ZZ0 = zeros(size(RR0));
pix4 = [RR0(:)'; CC0(:)'; ZZ0(:)'; ones(1,numel(RR0))];
XYZw = AF * pix4;               % 4xN, patient/world coords

ijk = AM_r \ XYZw;              % 4xN, index coords (0-based)
rowq = reshape(ijk(1,:) + 1, nR, nC);  % MATLAB 1-based
colq = reshape(ijk(2,:) + 1, nR, nC);
slcq = reshape(ijk(3,:) + 1, nR, nC);

%% Sample MPRAGE on the FASE grid
matched_slice = interpn( ...
    1:size(volM_r,1), 1:size(volM_r,2), 1:size(volM_r,3), double(volM_r), ...
    rowq, colq, slcq, 'linear', 0);

inside = rowq>=1 & rowq<=size(volM_r,1) & colq>=1 & colq<=size(volM_r,2) & slcq>=1 & slcq<=size(volM_r,3);
fprintf('Inside fraction: %.3f\n', nnz(inside)/numel(inside));

%% QA figure
fase_show = volF(:,:,1);
figure('Color','w','Position',[100,100,900,400]);
subplot(1,2,1); imshow(mat2gray(fase_show));     title('FASE 2D (original grid)');
subplot(1,2,2); imshow(mat2gray(matched_slice)); title('MPRAGE → FASE geometry');
exportgraphics(gcf,'Slice_Match.png'); close;
mkdir('Image4Segmentation')
cd 'Image4Segmentation'
dicomwrite(matched_slice,'IMG-0001-0001.dcm',infoF(1))
cd ..


end