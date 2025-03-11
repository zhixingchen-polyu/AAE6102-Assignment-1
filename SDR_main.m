%Purpose:
%   Main function of the software-defined radio (SDR) receiver platform
%
%--------------------------------------------------------------------------
%                           SoftXXXGPS v1.0
% 
% Copyright (C) X X  
% Written by X X


% 
clear; 
format long g;
addpath geo             %  
addpath acqtckpos       % Acquisition, tracking, and postiong calculation functions


%% Parameter initialization 
[file, signal, acq, track, solu, cmn] = initParameters();

 
%% Acquisition 
if ~exist(['Acquired_',file.fileName,'_',num2str(file.skip),'.mat'])
    Acquired = acquisition_hs(file,signal,acq); %
    save(['Acquired_',file.fileName,'_',num2str(file.skip)],'Acquired');    
else
    load(['Acquired_',file.fileName,'_',num2str(file.skip),'.mat']);
end 
fprintf('Acquisition Completed. \n\n');
 
% return


%% Do conventional signal tracking and obtain satellites ephemeris
fprintf('Tracking ... \n\n');
if ~exist(['eph_',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat'])
    % tracking using conventional DLL and PLL
    if ~exist(['TckResult_Eph',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']) %
        [TckResultCT, CN0_Eph] =  trackingCT(file,signal,track,Acquired); 
        TckResult_Eph = TckResultCT;
        save(['TckResult_Eph',file.fileName,'_',num2str(track.msToProcessCT/1000)], 'TckResult_Eph','CN0_Eph'); 
        %% Plotting the correlation outputs and saving figures
        for svIndex = 1:length(Acquired.sv)  % Loop through all acquired satellites
            if exist('TckResultCT', 'var')  % Ensure TckResultCT variable exists
                figure;  % Create a new figure
                plot(TckResultCT(Acquired.sv(svIndex)).P_i, 'g', 'DisplayName', 'P_i'); hold on;
                plot(TckResultCT(Acquired.sv(svIndex)).P_q, 'b', 'DisplayName', 'P_q');
                plot(TckResultCT(Acquired.sv(svIndex)).E_i, 'r', 'DisplayName', 'E_i');
                plot(TckResultCT(Acquired.sv(svIndex)).L_i, 'k', 'DisplayName', 'L_i');
                title(['Correlation Outputs for Satellite ', num2str(Acquired.sv(svIndex))]);
                xlabel('Sample Index');
                ylabel('Correlation Value');
                legend show;
                hold off;

                % Save the figure
                saveFileName = sprintf('%s_Satellite_%d.png', file.fileName, Acquired.sv(svIndex));
                saveas(gcf, saveFileName);  % Save the current figure as a PNG file
            else
                error('TckResultCT does not exist.'); % Handle case where it does not exist
            end
        end
    else   
        load(['TckResult_Eph',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']);
    end 
    
    % navigaion data decode
    fprintf('Navigation data decoding ... \n\n');
    [eph, ~, sbf] = naviDecode_updated(Acquired, TckResult_Eph);
    save(['eph_',file.fileName,'_',num2str(track.msToProcessCT/1000)], 'eph');
    save(['sbf_',file.fileName,'_',num2str(track.msToProcessCT/1000)], 'sbf');
%     save(['TckRstct_',file.fileName,'_',num2str(track.msToProcessCT/1000)], 'TckResultCT'); % Track results are revised in function naviDecode for 20 ms T_coh
else
    load(['eph_',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']);
    load(['sbf_',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']);
    load(['TckResult_Eph',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']);
end 

  
%% Find satellites that can be used to calculate user position
posSV  = findPosSV(file,Acquired,eph);
 
%% Do positiong in conventional or vector tracking mode
cnslxyz = llh2xyz(solu.iniPos); % initial position in ECEF coordinate
 
if cmn.vtEnable == 1    
    fprintf('Positioning (VTL) ... \n\n');
  
    % load data to initilize VT
    load(['nAcquired_',file.fileName,'_',num2str(file.skip),'.mat']); % load acquired satellites that can be used to calculate position  
    Acquired = nAcquired;  
    
    load(['eph_',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']); % load eph
    load(['sbf_',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']); % 
    
    load(['tckRstCT_1ms_',file.fileName,'.mat']);%,'_Grid'
    load(['navSolCT_1ms_',file.fileName,'.mat']); 
     
    [TckResultVT, navSolutionsVT] = ...
                  trackingVT_POS_updated(file,signal,track,cmn,solu,Acquired,cnslxyz,eph,sbf,TckResult_Eph, TckResultCT_pos,navSolutionsCT);  
else 
    load(['nAcquired_',file.fileName,'_',num2str(file.skip),'.mat']); % load acquired satellites that can be used to calculate position  
    Acquired = nAcquired;
    
    [TckResultCT_pos, navSolutionsCT] = ...
           trackingCT_POS_updated(file,signal,track,cmn,Acquired,TckResult_Eph, cnslxyz,eph,sbf,solu); %trackingCT_POS_multiCorr_1ms
                 
end 

fprintf('Tracking and Positioing Completed.\n\n');


 

