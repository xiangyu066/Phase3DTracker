%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: SwimmerSpatialTracker.m
% - Author: XYZ
% - Created date: April 9, 2020
% - Modified date: March 9, 2022
% - Notes:
%       1.) This version have to manually choose target, and the target
%       have to roughly put in the middle.
% - Next modified:
%       1.) Automatically chosen target
% - Version: 3.0
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, warning('off')
disp('Running...')

%% Define units
global um px sec msec
um = 1;
px = 1;
sec = 1;
msec = 1E-3 *(sec);

%% Define parameters of imaging system
% define imaging system parameters 
dz = 0.05*(um);                                                             % the axial depth between layers
pixelsize = 6.5*(um);
Obj_Mag = 40;                                                               % the magnification of objective
dt = 30*(msec);                                                             % the timestep between images

% image processing parameters
speed_upper = 120*(um/sec);
sz = -45*(um);                                                              % the lower working depth
ez = 45*(um);                                                               % the upper working depth

%
inputdir = 'E:\20220323';                                          % the directory path of the time-series movie
outputdir = 'E:\20220323';                                         % save tracking results

%% Preallocating 1vavriables and functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nFile = 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

listing = dir([inputdir,'\*.tif']);
inputfile = [inputdir,'\',listing(nFile).name];
nFrames = length(imfinfo(inputfile));                                       % the number of stacked images

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sFrame = 744                                                            % the starting tracking frame
eFrame = nFrames                                                            % the ending tracking frame
eFrame2 = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load depth library
load('E:\20210610\Lib\Library.mat')
load('E:\20210610\Lib\Library_Pos.mat')

% constrain library depth range
Library(:,:,Library_Pos<sz) = [];
Library_Pos(Library_Pos<sz) = [];
Library(:,:,Library_Pos>ez) = [];
Library_Pos(Library_Pos>ez) = [];
Lib_Square = size(Library,1);
nLayers = size(Library,3);

% normalized library
Library = (Library-mean(mean(Library,2),1))./std(Library,0,[1,2]);

% Preallocating vavriables
eff_pixelsize = pixelsize/Obj_Mag;
search_Box_xy = Lib_Square+2*round(speed_upper*dt/eff_pixelsize);
search_Box_z = 2*round(speed_upper*dt/dz);

% padding library
Lib_padding = zeros([search_Box_xy,search_Box_xy,nLayers]);
Lib_padding(1:Lib_Square,1:Lib_Square,:) = Library;

%
maxCrr_nLayer = 0;

% Make database folder
[~,msg] = mkdir([outputdir,'\Analyzed']);

%% Human-made define target
nSeeds = 0;
listSeeds = [];

origina_ = double(imread(inputfile,sFrame));
origina_(1:2,:) = 0;
origina_(end-1:end,:) = 0;
origina = zeros(size(origina_)+search_Box_xy);
origina(search_Box_xy/2+1:search_Box_xy/2+size(origina_,1),search_Box_xy/2+1:search_Box_xy/2+size(origina_,2)) = origina_;
figure(20),clf(gcf),imshow(origina,[mean2(origina_)-4*std2(origina_),mean2(origina_)+4*std2(origina_)]),
title([num2str(nFile),' / ',num2str(length(listing))])

nCounts = str2double(cell2mat(inputdlg('Enter number:')));
if (nCounts>0)
    for nCount = 1:nCounts
        nSeeds = nSeeds+1;
        h = drawrectangle('Label',['Target',num2str(nSeeds)],...
            'DrawingArea',[1,1,flip(size(origina))],...
            'Position',[1,1,search_Box_xy,search_Box_xy]);
        wait(h);
        position = h.Position;
        sx = round(position(1));
        sy = round(position(2));
        ex = sx+search_Box_xy-1;
        ey = sy+search_Box_xy-1;
        listSeeds(nSeeds,:) = [nFile,nCount,sx,sy,ex,ey];
    end
end

%% Tracking
figure(10), set(gcf,'WindowStyle','docked')
figure(30), set(gcf,'WindowStyle','docked')
pause(1)
for nSeed = 1:nSeeds
    % next saved file index
    nCount = 0;
    TF = true;
    while (TF)
        nCount = nCount+1;
        outputname = [strrep(listing(listSeeds(nSeed,1)).name,'.tif',''),'-',num2str(nCount)];
        TF = isfile([outputdir,'\Analyzed\',outputname,'.mat']);
    end
    
    % record a video of the tracking process
    outputname = [strrep(listing(listSeeds(nSeed,1)).name,'.tif',''),'-',num2str(nCount)]
    vSession = VideoWriter([outputdir,'\Analyzed\',outputname,'.avi']);
    vSession.FrameRate = 1/dt;
    open(vSession)
    
    Pos = NaN(nFrames,3);
    clf(10), tic
    % forward tracking
    for nFrame = sFrame:eFrame
        % read image
        inputfile = [inputdir,'\',listing(listSeeds(nSeed,1)).name];
        origina_ = double(imread(inputfile,nFrame));
        origina_(1:2,:) = 0;
        origina_(end-1:end,:) = 0;
        origina = zeros(size(origina_)+search_Box_xy);
        origina(search_Box_xy/2+1:search_Box_xy/2+size(origina_,1),search_Box_xy/2+1:search_Box_xy/2+size(origina_,2)) = origina_;

        % label tracking seeds
        if (nFrame==sFrame)
            sx = listSeeds(nSeed,3);
            sy = listSeeds(nSeed,4);
            ex = listSeeds(nSeed,5);
            ey = listSeeds(nSeed,6);
            sLayer = 1;
            eLayer = nLayers;
        else
            if ~isempty(maxCrr_nLayer)
                sx = sx+shiftX+Lib_Square/2-search_Box_xy/2;
                sy = sy+shiftY+Lib_Square/2-search_Box_xy/2;
                ex = sx+search_Box_xy-1;
                ey = sy+search_Box_xy-1;
            end
            sLayer = old_nLayer-search_Box_z/2;
            eLayer = old_nLayer+search_Box_z/2;
        end
        
        % determinate whether out of boundary
        checksum_xy = (sx>=10)+(sy>=10)+(ex<=size(origina,2)-10)+(ey<=size(origina,1)-10);
        checksum_z = (sLayer>=1)+(eLayer<=nLayers);
        if (checksum_xy<4) || (checksum_z<2)
            disp('Out of bound.');
            break;
        end
        
        % normalized
        neworigina = origina(sy:ey,sx:ex);
        neworigina(neworigina~=0) = (neworigina(neworigina~=0)-mean2(neworigina(neworigina~=0)))./std2(neworigina(neworigina~=0));
        
        % using FFT to boost calculating correlation efficiency
        Crrs = zeros(nLayers,4);
        for nLayer = sLayer:eLayer
            Crr = abs(ifft2(fft2(neworigina).*conj(fft2(Lib_padding(:,:,nLayer)))));
            [shiftY,shiftX] = find(Crr==max(Crr(:)));
            Crrs(nLayer,:) = [shiftX,shiftY,nLayer,max(Crr(:))];
        end
        
        % filt out of boundary
        if (nFrame>1)
            Crrs((Crrs(:,3)==0),:) = '';
            Crrs((Crrs(:,1)>size(neworigina,2)-Lib_Square+1),:) = '';
            Crrs((Crrs(:,2)>size(neworigina,1)-Lib_Square+1),:) = '';
        end
        
        % find maximal correlation coefficient
        maxCrr_nLayer = find(Crrs(:,4)==max(Crrs(:,4)));
        if isempty(maxCrr_nLayer)
            shiftX = 0;
            shiftY = 0;
        else
            shiftX = Crrs(maxCrr_nLayer,1);
            shiftY = Crrs(maxCrr_nLayer,2);
            old_shiftX = shiftX;
            old_shiftY = shiftY;
            old_nLayer = Crrs(maxCrr_nLayer,3);
        end
        
        % write databse
        if isempty(maxCrr_nLayer)
            Pos(nFrame,:) = [[sx+shiftX,sy+shiftY]*eff_pixelsize,Pos(nFrame-1,3)];
        else
            Pos(nFrame,:) = [[sx+shiftX,sy+shiftY]*eff_pixelsize,Library_Pos(Crrs(maxCrr_nLayer,3))];
        end
        
        % real-time monitoring
        figure(10),
        subplot 221, cla(gca), imshow(neworigina,[]),title(['nFrame = ',num2str(nFrame),' (',num2str((nFrame-1)*dt),' sec)'])
        hold on,rectangle('position',[old_shiftX,old_shiftY,Lib_Square,Lib_Square],'edgecolor','r')
        subplot 222, cla(gca), plot3(Pos(sFrame,1),Pos(sFrame,2),Pos(sFrame,3),'>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
        hold on,plot3(Pos(:,1),Pos(:,2),Pos(:,3),'k-')
        hold on,plot3(Pos(nFrame,1),Pos(nFrame,2),Pos(nFrame,3),'s','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
        xlabel('X [\mum]','fontweight','bold'),ylabel('Y [\mum]','fontweight','bold'),zlabel('Z [\mum]','fontweight','bold')
        grid on
        subplot 223, cla(gca), imshow(neworigina(old_shiftY:old_shiftY+Lib_Square-1,old_shiftX:old_shiftX+Lib_Square-1),[]),title('Dock target')
        hold on,line([Lib_Square/2+1,Lib_Square/2+1],[1,Lib_Square],'color','r')
        line([1,Lib_Square],[Lib_Square/2+1,Lib_Square/2+1],'color','r')
        subplot 224, cla(gca), imshow(Library(:,:,old_nLayer),[]),title(['nLayer = ',num2str(old_nLayer),' (',num2str(Library_Pos(old_nLayer)),' \mum)'])
        hold on,line([Lib_Square/2+1,Lib_Square/2+1],[1,Lib_Square],'color','r')
        line([1,Lib_Square],[Lib_Square/2+1,Lib_Square/2+1],'color','r')
        drawnow
        
        % save tracking process
        frame = getframe(gcf);
        writeVideo(vSession,frame);
    end
    
    % backward tracking
    for nFrame = sFrame:-1:eFrame2
        % read image
        inputfile = [inputdir,'\',listing(listSeeds(nSeed,1)).name];
        origina_ = double(imread(inputfile,nFrame));
        origina_(1:2,:) = 0;
        origina_(end-1:end,:) = 0;
        origina = zeros(size(origina_)+search_Box_xy);
        origina(search_Box_xy/2+1:search_Box_xy/2+size(origina_,1),search_Box_xy/2+1:search_Box_xy/2+size(origina_,2)) = origina_;

        % label tracking seeds
        if (nFrame==sFrame)
            sx = listSeeds(nSeed,3);
            sy = listSeeds(nSeed,4);
            ex = listSeeds(nSeed,5);
            ey = listSeeds(nSeed,6);
            sLayer = 1;
            eLayer = nLayers;
        else
            if ~isempty(maxCrr_nLayer)
                sx = sx+shiftX+Lib_Square/2-search_Box_xy/2;
                sy = sy+shiftY+Lib_Square/2-search_Box_xy/2;
                ex = sx+search_Box_xy-1;
                ey = sy+search_Box_xy-1;
            end
            sLayer = old_nLayer-search_Box_z/2;
            eLayer = old_nLayer+search_Box_z/2;
        end
        
        % determinate whether out of boundary
        checksum_xy = (sx>=10)+(sy>=10)+(ex<=size(origina,2)-10)+(ey<=size(origina,1)-10);
        checksum_z = (sLayer>=1)+(eLayer<=nLayers);
        if (checksum_xy<4) || (checksum_z<2)
            disp('Out of bound.');
            break;
        end
        
        % normalized
        neworigina = origina(sy:ey,sx:ex);
        neworigina(neworigina~=0) = (neworigina(neworigina~=0)-mean2(neworigina(neworigina~=0)))./std2(neworigina(neworigina~=0));
        
        % using FFT to boost calculating correlation efficiency
        Crrs = zeros(nLayers,4);
        for nLayer = sLayer:eLayer
            Crr = abs(ifft2(fft2(neworigina).*conj(fft2(Lib_padding(:,:,nLayer)))));
            [shiftY,shiftX] = find(Crr==max(Crr(:)));
            Crrs(nLayer,:) = [shiftX,shiftY,nLayer,max(Crr(:))];
        end
        
        % filt out of boundary
        if (nFrame>1)
            Crrs((Crrs(:,3)==0),:) = '';
            Crrs((Crrs(:,1)>size(neworigina,2)-Lib_Square+1),:) = '';
            Crrs((Crrs(:,2)>size(neworigina,1)-Lib_Square+1),:) = '';
        end
        
        % find maximal correlation coefficient
        maxCrr_nLayer = find(Crrs(:,4)==max(Crrs(:,4)));
        if isempty(maxCrr_nLayer)
            shiftX = 0;
            shiftY = 0;
        else
            shiftX = Crrs(maxCrr_nLayer,1);
            shiftY = Crrs(maxCrr_nLayer,2);
            old_shiftX = shiftX;
            old_shiftY = shiftY;
            old_nLayer = Crrs(maxCrr_nLayer,3);
        end
        
        % write databse
        if isempty(maxCrr_nLayer)
            Pos(nFrame,:) = [[sx+shiftX,sy+shiftY]*eff_pixelsize,Pos(nFrame-1,3)];
        else
            Pos(nFrame,:) = [[sx+shiftX,sy+shiftY]*eff_pixelsize,Library_Pos(Crrs(maxCrr_nLayer,3))];
        end
        
        % real-time monitoring
        figure(10),
        subplot 221, cla(gca), imshow(neworigina,[]),title(['nFrame = ',num2str(nFrame),' (',num2str((nFrame-1)*dt),' sec)'])
        hold on,rectangle('position',[old_shiftX,old_shiftY,Lib_Square,Lib_Square],'edgecolor','r')
        subplot 222, cla(gca), plot3(Pos(sFrame,1),Pos(sFrame,2),Pos(sFrame,3),'>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
        hold on,plot3(Pos(:,1),Pos(:,2),Pos(:,3),'k-')
        hold on,plot3(Pos(eFrame,1),Pos(eFrame,2),Pos(eFrame,3),'s','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
        hold on,plot3(Pos(nFrame,1),Pos(nFrame,2),Pos(nFrame,3),'o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
        xlabel('X [\mum]','fontweight','bold'),ylabel('Y [\mum]','fontweight','bold'),zlabel('Z [\mum]','fontweight','bold')
        grid on
        subplot 223, cla(gca), imshow(neworigina(old_shiftY:old_shiftY+Lib_Square-1,old_shiftX:old_shiftX+Lib_Square-1),[]),title('Dock target')
        hold on,line([Lib_Square/2+1,Lib_Square/2+1],[1,Lib_Square],'color','r')
        line([1,Lib_Square],[Lib_Square/2+1,Lib_Square/2+1],'color','r')
        subplot 224, cla(gca), imshow(Library(:,:,old_nLayer),[]),title(['nLayer = ',num2str(old_nLayer),' (',num2str(Library_Pos(old_nLayer)),' \mum)'])
        hold on,line([Lib_Square/2+1,Lib_Square/2+1],[1,Lib_Square],'color','r')
        line([1,Lib_Square],[Lib_Square/2+1,Lib_Square/2+1],'color','r')
        drawnow
        
        % save tracking process
        frame = getframe(gcf);
        writeVideo(vSession,frame);
    end
    
    close(vSession),toc
    
    % plot swimmer's trajectory
    figure(30),
    subplot 221, cla(gca), plot3(Pos(sFrame,1),Pos(sFrame,2),Pos(sFrame,3),'>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on,plot3(Pos(:,1),Pos(:,2),Pos(:,3),'k')
    hold on,plot3(Pos(eFrame,1),Pos(eFrame,2),Pos(eFrame,3),'s','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    xlabel('X [\mum]','fontweight','bold'),ylabel('Y [\mum]','fontweight','bold'),zlabel('Z [\mum]','fontweight','bold')
    grid on,axis equal
    subplot 222, cla(gca), plot(Pos(sFrame,1),Pos(sFrame,2),'>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on,plot(Pos(nFrame,1),Pos(nFrame,2),'s','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on,plot(Pos(:,1),Pos(:,2),'k')
    xlabel('X [\mum]','fontweight','bold'), ylabel('Y [\mum]','fontweight','bold')
    legend({'Start','End'})
    grid on, axis equal
    subplot 223, cla(gca), plot(Pos(sFrame,1),Pos(sFrame,3),'>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on,plot(Pos(nFrame,1),Pos(nFrame,3), 's','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on,plot(Pos(:,1),Pos(:,3),'k')
    xlabel('X [\mum]','fontweight','bold'),ylabel('Z [\mum]','fontweight','bold')
    legend({'Start','End'})
    grid on,axis equal
    subplot 224, cla(gca), plot(Pos(sFrame,2),Pos(sFrame,3),'>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on,plot(Pos(nFrame,2),Pos(nFrame,3),'s','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on,plot(Pos(:,2),Pos(:,3),'k')
    xlabel('Y [\mum]','fontweight','bold'),ylabel('Z [\mum]','fontweight','bold')
    legend({'Start','End'})
    grid on, axis equal
    frame = getframe(gcf);
    imwrite(frame.cdata,[outputdir,'\Analyzed\',outputname,'.png'])
    
    % save tracking data
    data.unit = {'um','sec'};
    data.dt = dt;
    data.Pos = Pos;
    save([outputdir,'\Analyzed\',outputname,'.mat'],'data')
end

%%
disp('Done.')
