echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: SwimmerSpatialTracker_LibBuilder.m
% - Author: XYZ
% - Created date: April 6, 2020
% - Modified date: June 17, 2021
% - Notes:
%       1.) This code is limited in "ONE frame only contains ONE cell".
% - Next modified:
%       1.) Adaptive threshold
% - Version: 2.2
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
close all, clear all, warning('off')
disp('Running...')

%%
global um px
um = 1;
px = 1;

%% Define image infomation
inputdir = 'Z:\PhD in NCU\Data\Experiments\Beads\PS beads\3D tracking\20210316\1um PS\Zstep_lib\1umPS_dz0.2um';
dz = 0.2 *(um);
Lib_Square = 200 *(px);                                                     % even number
pixelsize = 6.5 *(um);
M = 40;

%% Preallocating vavriables and functions
% extract information about graphics file
listing = dir([inputdir, '\*.tif']);
nFiles = length(listing);
info = imfinfo([inputdir, '\', listing(1).name]);
width = info(1).Width *(px);
height = info(1).Height *(px);
nLayers = length(info);
eff_pixelsize = pixelsize/M;

% Preallocating vavriables
Proj_Images = zeros(height,width,nFiles);
nCells = 0;                                                                 % count for availible cells

% Make database folder
% [~,msg] = mkdir('Lib0');
% if ~isempty(msg)
%     error('Please remove the existing "Lib" folder.')
% end

%% Do maximal projection for all chosen cells
disp('Step (1/5): Do maxial projection...'), tic

for nFile = 1:nFiles
    neworigina = zeros(height,width);
    disp(['Processing...(', num2str(nFile), '/', num2str(nFiles), ')']);
    for nLayer = 1:nLayers
        origina = double(imread([inputdir, '\', listing(nFile).name], nLayer));
        neworigina(origina>neworigina) = origina(origina>neworigina);
    end
    Proj_Images(:,:, nFile) = neworigina;
    
%     if (nFile == 1)
%         imwrite(uint16(Proj_Images(:,:, nFile)), [pwd, '\Lib\', 'Proj_Images.tif'], 'Compression', 'none')
%     else
%         imwrite(uint16(Proj_Images(:,:, nFile)), [pwd, '\Lib\', 'Proj_Images.tif'], 'WriteMode', 'append', 'Compression', 'none')
%     end
end
toc

%% Put the maximal intensity into the center of the library database
disp('Step (2/5): Analyze the centroid from maximal projection...'), tic

CMxys = zeros(nFiles,2);
idxs = [];
figure(1), set(gcf,'WindowStyle', 'docked')
for nFile = 1:nFiles
    disp(['Processing...(', num2str(nFile), '/', num2str(nFiles), ')']);
    
    cla(gca), imshow(Proj_Images(:,:, nFile),[]), title([num2str(nFile), '/', num2str(nFiles)])
    
    thresh = 0.9*max(max(Proj_Images(height/2-32:height/2+32,width/2-32:width/2+32, nFile)));
    mask = (Proj_Images(:,:,nFile) >= thresh);
    mask = bwareaopen(mask,4,8);
    s = regionprops(mask, 'centroid');
    CMxy = cat(1, s.Centroid);
    hold on, plot(CMxy(:,1), CMxy(:,2), 'gs')
    
    % filt out near image boundary
    checksum = sum((abs(CMxy-[width, height])>=[Lib_Square/2, Lib_Square/2]), 2)...
        + sum((abs(CMxy-[1, 1])>=[Lib_Square/2,Lib_Square/2]), 2);
    
    % record centroids
    if (checksum==4)
        nCells = nCells +1;
        CMxys(nCells,:) = round(CMxy);
        
        hold on, plot(CMxy(:,1),CMxy(:,2),'r+')
        frame = getframe(gcf);
%         imwrite(frame.cdata, [pwd,'\Lib\','CM_', strrep(listing(nFile).name,'.tif',''), '.png'])
    else
        idxs = [idxs,nFile];
    end
end
listing(idxs) = '';
disp(['There are ',num2str(nCells), ' available cells and filt out ', num2str(nFiles-nCells), ' cells.']), toc

%% Shift the maxial axial intensity into zero
disp('Step (3/5): Align layers in maxial axial intensity...'), tic

alignLayers = zeros(nLayers, nCells);
CentralVals = zeros(nLayers, nCells);
for nCell = 1:nCells
    disp(['Processing...(', num2str(nCell), '/', num2str(nCells), ')']);
    
    CMxy = CMxys(nCell,:);
    for nLayer = 1:nLayers
        origina = double(imread([inputdir,'\',listing(nCell).name], nLayer));
        CentralVals(nLayer, nCell) = origina(sub2ind([height, width],CMxy(2),CMxy(1)));
    end
    
    % align with the maximal central intensity
    zeros_Layer = find(CentralVals(:, nCell) == max(CentralVals(:, nCell)));
    alignLayers(:, nCell) = (1:nLayers)' - zeros_Layer(end);
end
toc
% [rows,~] = find(alignLayers==0);
% cols = find(abs(rows-round(nLayers/2))>200);
% nCells = nCells-cols;
% alignLayers(:,cols) = '';
% listing(cols) = '';
% disp(['There are ',num2str(nCells), ' available cells and filt out ', num2str(cols), ' cells.'])

%% Build a library
disp('Step (4/5): Build a library...'), tic

upper_Layer = max(alignLayers(:));
lower_Layer = min(alignLayers(:));
Lib_nCounts = zeros(upper_Layer-lower_Layer+1,1); % count for contributing cell number 
Library = zeros(Lib_Square, Lib_Square, upper_Layer-lower_Layer+1);
alignCentralVals = zeros(upper_Layer-lower_Layer+1, nCells);
Library_Pos = -dz*(lower_Layer:upper_Layer)';

for nCell = 1:nCells
    disp(['Processing...(', num2str(nCell), '/', num2str(nCells), ')']);
    
    CMxy = CMxys(nCell,:);
    for nLayer = 1:nLayers
        origina = double(imread([inputdir,'\',listing(nCell).name], nLayer));
        
        sx = CMxy(1)-Lib_Square/2;
        sy = CMxy(2)-Lib_Square/2;
        ex = CMxy(1)+Lib_Square/2-1;
        ey = CMxy(2)+Lib_Square/2-1;
        neworigina = origina(sy:ey,sx:ex);
        
        nLayer_idx = alignLayers(nLayer,nCell) - lower_Layer + 1;
        Library(:, :, nLayer_idx) = Library(:, :, nLayer_idx) + neworigina;
        alignCentralVals(nLayer_idx,nCell) = CentralVals(nLayer,nCell);
        Lib_nCounts(nLayer_idx) = Lib_nCounts(nLayer_idx)+1;
    end
end
Library = Library/nCells;
toc

% using rotation-averaged to remove non-uniform in xy plane because symm.
% Library = (imrotate(Library,90) + imrotate(Library,180) + imrotate(Library,270) + imrotate(Library,360))/4;

% remove layers if contributing cell number is less than expected number
Library(:,:,Lib_nCounts<nCells) = '';
alignCentralVals(Lib_nCounts<nCells,:) = '';
Library_Pos(Lib_nCounts<nCells) = '';

% save database
disp('Step (5/5): Save database...'), tic

for nLayer = 1:size(Library,3)
    if (nLayer == 1)
        imwrite(uint16(Library(:,:,nLayer)),[pwd,'\Lib\','Library.tif'],'Compression','none')
    else
        imwrite(uint16(Library(:,:,nLayer)),[pwd,'\Lib\','Library.tif'],'WriteMode','append','Compression','none')
    end
end
save([pwd,'\Lib\Library.mat'], 'Library')
save([pwd,'\Lib\Library_Pos.mat'], 'Library_Pos')
toc

%% show results
y = ((1:Lib_Square)-(Lib_Square/2+1))*eff_pixelsize;
ref_x0 = Lib_Square/2+1;
ref_y0 = Lib_Square/2+1;
for nLayer = 1:size(Library,3)
    ref_x(:,nLayer) = Library(:,ref_x0,nLayer);
    ref_y(nLayer,:) = Library(ref_y0,:,nLayer);
end
figure(2), set(gcf,'WindowState','maximized')
subplot 211, imagesc(Library_Pos,y,log2(ref_x)), colormap('jet'), colorbar
xlabel('Z-position [\mum]','fontsize',24,'fontweight','bold')
ylabel('X-position [\mum]','fontsize',24,'fontweight','bold')
set(gca,'fontsize',16)
xlim([-50,50])
subplot 212, imagesc(Library_Pos,y,log2(ref_y.')), colormap('jet'), colorbar
xlabel('Z-position [\mum]','fontsize',24,'fontweight','bold')
ylabel('Y-position [\mum]','fontsize',24,'fontweight','bold')
set(gca,'fontsize',16)
xlim([-50,50])
frame = getframe(gcf);
% imwrite(frame.cdata, [pwd,'\Lib\','The response curve.png'])

figure(3), set(gcf,'WindowState','maximized')
plot(repmat(Library_Pos,[nCells,1]),alignCentralVals(:),'c.')
hold on, plot(Library_Pos, mean(alignCentralVals,2),'r')
legend({['nCells = ',num2str(nCells)],'Averaged'})
xlabel('Z [\mum]','fontweight','bold','fontsize',24), 
ylabel('The central intensity [a.u.]','fontweight','bold','fontsize',24)
set(gca,'fontsize',16)
xlim([-50,50])
grid on
frame = getframe(gcf);
% imwrite(frame.cdata, [pwd,'\Lib\','The axial profile of library.png'])

%%
disp('Done.')

