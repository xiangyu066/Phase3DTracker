%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Title: SwimmerSpatialTracker_MotionAnalyzer.m
% - Author: XYZ
% - Created date: April 9, 2020
% - Modified date: March 19, 2022
% - Notes:
%       1.) 
% - Next modified:
%       1.) 
% - Version: 2.9
% - Environments: Win10 (64-bit) / MATLAB 2019a (64-bit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, warning('off')
disp('Running...'), tic

%% Define units
global um msec sec
um = 1;
sec = 1;
msec = 1E-3 *(sec);

%% 
inputdir = 'G:\我的雲端硬碟\Data\ConstantEnv_TMN\2mM Mg2+\300mM Na+\20210929-2'
order = 2;
framelen = 7;
framelen_align = 3;
bundle_win = floor(framelen/2);
vMax = 100*(um/sec);
isSaveFig = false

%% load file
listing = dir([inputdir,'\*.mat']);
nFiles = length(listing)
outputdir = ''

%%
turn_angles = [];
turn_speeds = [];
for nFile = 9%1:nFiles
    % initial variables
    inputfile = [inputdir,'\',listing(nFile).name]
    load(inputfile);
    
    % extract data
    dt = (data.dt) *(sec);
    Pos = (data.Pos) *(um);
    
    % remove nan-term
    checknan = isnan(Pos(:,1)) & isnan(Pos(:,2)) & isnan(Pos(:,3));
    Pos(checknan,:) = [];
    
    % Initialization
    time = (0:size(Pos,1)-1)*dt;
    swimming_speeds = NaN(size(time));
    align_T = NaN(length(time),1);
    
    % Analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define first point as reference point
    x = Pos(:,1)-Pos(1,1);
    y = Pos(:,2)-Pos(1,2);
    z = Pos(:,3)-Pos(1,3);
    
    % Savitzky-Golay filtering
    new_x = sgolayfilt(x,order,framelen);                                   % output: new_x, new_y, new_z
    new_y = sgolayfilt(y,order,framelen);
    new_z = sgolayfilt(z,order,framelen);
    
    % calculate swimming speed
    for nPt = 2:size(Pos,1)-1
        dr = [new_x(nPt+1)-new_x(nPt-1),new_y(nPt+1)-new_y(nPt-1),new_z(nPt+1)-new_z(nPt-1)];
        swimming_speed = sqrt(sum(dr.^2))/(2*dt);
        swimming_speeds(nPt) = swimming_speed;                              % output: swimming_speeds
    end
    
    % calculate frenet
    [T,N,B,k,t] = frenet(new_x,new_y,new_z);                                % output: T, N, B, k, t
    for i = 1+bundle_win:size(T,1)-bundle_win
        align_T(i,1) = dot(T(i,:), T(i+bundle_win,:));
    end
    
    % correct path alignment
    swim_dir_ = real(acosd(align_T));                                       % temperatory swimming direction  
    time_ = (1:length(swim_dir_))*dt;                                       % temperory time
    swim_dir = swim_dir_;
    pi_tf = (swim_dir_>90);                                                 % determine reverse direction
    zero_tf = (swim_dir_<90);                                               % determine same direction
    regime_pi = {};                                                         % label reverse regime
    regime_zero = {};                                                       % label same regime
    regime_pi_ = [];
    regime_zero_ = [];
    regime_pi_count = 0;
    regime_zero_count = 0;
    for i = 1:length(swim_dir_)
        if zero_tf(i) == 1
            regime_zero_ = [regime_zero_,i];
            if zero_tf(i+1)==0
                regime_zero_count = regime_zero_count+1;
                regime_zero{regime_zero_count} = regime_zero_;
                regime_zero_ = [];
            end
        end
        
        if pi_tf(i) == 1
            regime_pi_ = [regime_pi_,i];
            if pi_tf(i+1)==0
                regime_pi_count = regime_pi_count+1;
                regime_pi{regime_pi_count} = regime_pi_;
                regime_pi_ = [];
            end
        end
    end
    
    for i = 2:2:length(regime_pi)
        regime_pi_ = regime_pi{i};
        swim_dir(regime_pi_) = 180-swim_dir(regime_pi_);
    end
    
    for i = 2:2:length(regime_zero)
        regime_zero_ = regime_zero{i};
        swim_dir(regime_zero_) = 180-swim_dir(regime_zero_);
    end
    
    swim_dir = sgolayfilt(swim_dir,order,framelen_align);                   % smooth alignment signal; output: swim_dir
    
    % calculate angular speed
    angular_speeds = swim_dir_*(pi/180)/(bundle_win*dt);                    % output: angular_speeds
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure(10), set(1,'WindowState','maximized')
%     plot(time_,swim_dir_,'--',time_,swim_dir)
%     legend({'Before','After'})
%     set(gca,'fontsize',16)
%     xlabel('Elapsed time [sec]','fontsize',24,'fontweight','bold')
%     ylabel('Path alignment [\circ]','fontsize',24,'fontweight','bold')
%     grid on
    
    % Plot trajectory, speed, and path alignment %%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(20), set(20,'WindowState','maximized')
    subplot(2,3,[1,2,4,5]), cla(gca)
    surface([new_x,new_x]',[new_y,new_y]',[new_z,new_z]',[swimming_speeds;swimming_speeds],...
        'facecolor','none','edgecolor','interp','linewidth',2);
    hold on, plot3(new_x(1),new_y(1),new_z(1),'>','Color','r','MarkerSize',12,'MarkerFaceColor','#D9FFFF')
    hold on,plot3(new_x(end),new_y(end),new_z(end),'s','Color','b','MarkerSize',12,'MarkerFaceColor','#D9FFFF')
    set(gca,'fontsize',16)
    legend({'','Start','End'})
    xlabel('X [\mum]','fontsize',24,'fontweight','bold')
    ylabel('Y [\mum]','fontsize',24,'fontweight','bold')
    zlabel('Z [\mum]','fontsize',24,'fontweight','bold')
    colormap('jet'), c = colorbar; %caxis([0 vMax])
    c.Label.String = 'Speed [\mum \cdot sec^{-1}]';
    c.Label.FontSize = 24;
    c.Label.FontWeight = 'bold';
    daspect([1,1,1])
    grid on, view([-37.5,30])
    
    subplot(2,3,3), cla(gca)
    set(gca,'fontsize',16)
    yyaxis left, plot(time,swimming_speeds,'.-')
    ylim([0,vMax])
    xlabel('Elapsed time [sec]','fontsize',24,'fontweight','bold')
    ylabel('Speed [\mum \cdot sec^{-1}]','fontsize',24,'fontweight','bold')
    yyaxis right, plot(time,cosd(swim_dir))
    ylim([-1.05,1.05])
    xlim([0,max(time)])
    ylabel('Path alignment, Cos(\theta)','fontsize',24,'fontweight','bold')
    box on, grid on
    
    subplot(2,3,6), cla(gca)
    histogram(swimming_speeds,'binwidth',2.5*(um/sec),'Normalization','probability')
    set(gca,'fontsize',16)
    xlabel('Speed [\mum \cdot sec^{-1}]','fontsize',24,'fontweight','bold')
    xlim([0,vMax])
    xticks([0:20:vMax])
    ylabel('Probability','fontsize',24,'fontweight','bold')
    grid on
    
    if isSaveFig==true
        saveas(gcf,[outputdir,'\Speed\',strrep(listing(nFile).name,'.mat',''),'_Speed.png'])
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_slope=NaN(size(time));
    y_slope=NaN(size(time));
    z_slope=NaN(size(time));
    switch_timePts = zeros(size(time));
    x_switch_timePts = [];
    y_switch_timePts = [];
    z_switch_timePts = [];
    span = 2;%floor(framelen/2);
    for i = 1+span-1:length(new_x)-span+1
        if i>=1+span && i<=length(new_x)-span
            p_x = polyfit(time(i-span:i+span),new_x(i-span:i+span),1);
            p_y = polyfit(time(i-span:i+span),new_y(i-span:i+span),1);
            p_z = polyfit(time(i-span:i+span),new_z(i-span:i+span),1);
        else
            p_x = polyfit(time(i-span+1:i+span-1),new_x(i-span+1:i+span-1),1);
            p_y = polyfit(time(i-span+1:i+span-1),new_y(i-span+1:i+span-1),1);
            p_z = polyfit(time(i-span+1:i+span-1),new_z(i-span+1:i+span-1),1);
        end
        
        % check x whther turn
        if abs(p_x(1))>3*(um/sec)
            if p_x(1)>0
                x_slope(i)=1;
            else
                x_slope(i)=-1;
            end
        else
            x_slope(i)=0;
        end
        
        if i==1+span-1
            x_switch_state = x_slope(i);
        else
           if x_slope(i)~=0 
               if abs(x_slope(i)-x_switch_state)==2
                   x_switch_state = x_slope(i);
                   x_switch_timePts = [x_switch_timePts,i];
                   if i-span<1
                       switch_timePts(1:i) = 1;
                   else
                       switch_timePts(i-span:i) = 1;
                   end
                   
                   if i+span>length(time)
                       switch_timePts(i:end) = 1;
                   else
                       switch_timePts(i:i+span) = 1;
                   end
               else
                   x_switch_state = x_slope(i);
               end
           end
        end
        
        % check y whther turn
        if abs(p_y(1))>3*(um/sec)
            if p_y(1)>0
                y_slope(i)=1;
            else
                y_slope(i)=-1;
            end
        else
            y_slope(i)=0;
        end
        
        if i==1+span-1
            y_switch_state = y_slope(i);
        else
            if y_slope(i)~=0
                if abs(y_slope(i)-y_switch_state)==2
                    y_switch_state = y_slope(i);
                    y_switch_timePts = [y_switch_timePts,i];
                    if i-span<1
                        switch_timePts(1:i) = 1;
                    else
                        switch_timePts(i-span:i) = 1;
                    end
                    
                    if i+span>length(time)
                        switch_timePts(i:end) = 1;
                    else
                        switch_timePts(i:i+span) = 1;
                    end
                else
                    y_switch_state = y_slope(i);
                end
            end
        end
        
        % check z whther turn
        if abs(p_z(1))>5*(um/sec)
            if p_z(1)>0
                z_slope(i)=1;
            else
                z_slope(i)=-1;
            end
        else
            z_slope(i)=0;
        end
        
        if i==1+span-1
            z_switch_state = z_slope(i);
        else
            if z_slope(i)~=0 
                if abs(z_slope(i)-z_switch_state)==2
                    z_switch_state = z_slope(i);
                    z_switch_timePts = [z_switch_timePts,i];
                    if i-span<1
                        switch_timePts(1:i) = 1;
                    else
                        switch_timePts(i-span:i) = 1;
                    end
                    
                    if i+span>length(time)
                        switch_timePts(i:end) = 1;
                    else
                        switch_timePts(i:i+span) = 1;
                    end
                else
                    z_switch_state = z_slope(i);
                end
            end
        end
    end
    
    % collect turning angle
    [L,nSegs] = bwlabel(switch_timePts);
    for nSeg=1:nSegs
        if sum(angular_speeds(L==nSeg)>10)>0
            turn_angles=[turn_angles,max(swim_dir_(L==nSeg))];
            turn_speeds=[turn_speeds,max(angular_speeds(L==nSeg))];
        else
            switch_timePts(L==nSeg)=0;
        end
    end
    
    figure(60), set(60,'WindowState','maximized')
    subplot 231, cla(gca)
    plot(new_x,new_y)
    hold on, plot(new_x(1),new_y(1),'>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on, plot(new_x(end),new_y(end),'s','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on, plot(new_x(switch_timePts==1),new_y(switch_timePts==1),'r.')
    set(gca,'fontsize',16)
    xlabel('X [\mum]','fontsize',24,'fontweight','bold')
    ylabel('Y [\mum]','fontsize',24,'fontweight','bold')
    legend({'Trajectory','Start','Goal','Turning region'},'Location','best')
    daspect([1,1,1]), grid on
    
    subplot 232, cla(gca)
    plot(new_x,new_z)
    hold on, plot(new_x(1),new_z(1),'>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on, plot(new_x(end),new_z(end),'s','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on, plot(new_x(switch_timePts==1),new_z(switch_timePts==1),'r.')
    set(gca,'fontsize',16)
    xlabel('X [\mum]','fontsize',24,'fontweight','bold')
    ylabel('Z [\mum]','fontsize',24,'fontweight','bold')
    daspect([1,1,1]), grid on
    
    subplot 233, cla(gca)
    plot(new_y,new_z)
    hold on, plot(new_y(1),new_z(1),'>','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on, plot(new_y(end),new_z(end),'s','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    hold on, plot(new_y(switch_timePts==1),new_z(switch_timePts==1),'r.')
    set(gca,'fontsize',16)
    xlabel('Y [\mum]','fontsize',24,'fontweight','bold')
    ylabel('Z [\mum]','fontsize',24,'fontweight','bold')
    daspect([1,1,1]), grid on
    
    subplot 234, cla(gca)
    yyaxis left, plot(time,new_x,'.-',time(x_switch_timePts),new_x(x_switch_timePts),'mo',time(switch_timePts==1),new_x(switch_timePts==1),'gs'), 
    ylabel('X [um]','fontweight','bold','fontsize',16), 
    yyaxis right, plot(time,x_slope), ylim([-1.05,1.05]), grid on
    subplot 235, cla(gca) 
    yyaxis left, plot(time,new_y,'.-',time(y_switch_timePts),new_y(y_switch_timePts),'mo',time(switch_timePts==1),new_y(switch_timePts==1),'gs'), 
    ylabel('Y [um]','fontweight','bold','fontsize',16)
    yyaxis right, plot(time,y_slope), ylim([-1.05,1.05]), grid on
    subplot 236, cla(gca)
    yyaxis left, plot(time,new_z,'.-',time(z_switch_timePts),new_z(z_switch_timePts),'mo',time(switch_timePts==1),new_z(switch_timePts==1),'gs'), 
    xlabel('Time [sec]','fontweight','bold','fontsize',16), ylabel('Z [um]','fontweight','bold','fontsize',16), 
    yyaxis right, plot(time,z_slope), ylim([-1.05,1.05])
    ylabel('Slope','fontweight','bold'), grid on
    
    if isSaveFig==true
        saveas(gcf,[outputdir,'\Turning angle\',strrep(listing(nFile).name,'.mat',''),'_Turning.png'])
    end
    
   
    % Plot frenet vector for local trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     sidx = 1;
%     eidx = length(new_x);
%     
%     figure(70), set(70,'WindowState','maximized')
%     subplot(2,3,3), cla(gca)
%     surface([new_x(sidx:eidx),new_x(sidx:eidx)]',...
%         [new_y(sidx:eidx),new_y(sidx:eidx)]',...
%         [new_z(sidx:eidx),new_z(sidx:eidx)]',...
%         [swimming_speeds(sidx:eidx);swimming_speeds(sidx:eidx)],...
%         'facecolor','none','edgecolor','interp','linewidth',2);
%     set(gca,'fontsize',16)
%     xlabel('X [\mum]','fontsize',24,'fontweight','bold')
%     ylabel('Y [\mum]','fontsize',24,'fontweight','bold')
%     zlabel('Z [\mum]','fontsize',24,'fontweight','bold')
%     colormap('jet'), c= colorbar(); caxis([0 vMax])
%     c.Label.String = 'Speed [\mum \cdot sec^{-1}]';
%     c.Label.FontSize = 24;
%     c.Label.FontWeight = 'bold';
%     box on, grid on, view([-37.5,30])
%     
%     subplot(2,3,[1,2,4,5]), cla(gca)
%     plot3(new_x(sidx:eidx),new_y(sidx:eidx),new_z(sidx:eidx),'ko')
%     hold on,
%     quiver3(new_x(sidx:eidx),new_y(sidx:eidx),new_z(sidx:eidx),...
%         T(sidx:eidx,1),T(sidx:eidx,2),T(sidx:eidx,3),0,'color','r','linewidth',1)
%     quiver3(new_x(sidx:eidx),new_y(sidx:eidx),new_z(sidx:eidx),...
%         N(sidx:eidx,1),N(sidx:eidx,2),N(sidx:eidx,3),0,'color','g','linewidth',1)
%     quiver3(new_x(sidx:eidx),new_y(sidx:eidx),new_z(sidx:eidx),...
%         B(sidx:eidx,1),B(sidx:eidx,2),B(sidx:eidx,3),0,'color','b','linewidth',1)
%     set(gca,'fontsize',16)
%     xlabel('X [\mum]','fontsize',24,'fontweight','bold')
%     ylabel('Y [\mum]','fontsize',24,'fontweight','bold')
%     zlabel('Z [\mum]','fontsize',24,'fontweight','bold')
%     box on, grid on
%     
%     subplot(2,3,6), cla(gca)
%     yyaxis left, plot(time(sidx:eidx),swimming_speeds(sidx:eidx),'.-')
%     ylim([0,vMax])
%     set(gca,'fontsize',16)
%     xlabel('Elapsed time [sec]','fontsize',24,'fontweight','bold')
%     ylabel('Speed [\mum \cdot sec^{-1}]','fontsize',24,'fontweight','bold')
%     yyaxis right, plot(time_(sidx:eidx-1),cosd(swim_dir(sidx:eidx-1)))
%     ylim([-1.05,1.05]), grid on
%     ylabel('Path alignment, Cos(\theta)','fontsize',24,'fontweight','bold')
   
end

% plot histogram of turning angle
figure(80), h=histogram(turn_angles,'binwidth',20);
x=[];
y=[];
y=h.Values;
x_=h.BinEdges;
for i = 1:length(x_)-1
    x(i)=(x_(i)+x_(i+1))/2;
end
options = fitoptions('gauss2','Lower',[-Inf 30 0 -Inf 40 0],'Upper',[Inf 180 Inf Inf 180 Inf],'StartPoint',[1,180,40,1,70,20]);
f = fit(x.',y.','gauss2',options)
hold on, plot(f,x,y)
legend('Location','Best')
xlim([0,180]), xticks([0:45:180]), grid on
xlabel('Turning angle [deg]','fontweight','bold','fontsize',16)
ylabel('Count','fontweight','bold','fontsize',16)
set(gca,'fontsize',12)

figure(81), histogram(turn_speeds,'binwidth',4)
xlim([0,35])
grid on
xlabel('Turning speed [rad \cdot sec^{-1}]','fontweight','bold')
ylabel('Count','fontweight','bold')

%%
toc, disp('Done.')
