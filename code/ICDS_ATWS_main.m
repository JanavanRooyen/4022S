%% Image Contrast and Doppler Spread based Automatic Time Window Selection
% adjust parameters according to which data to save, and where to.

clear all
clc
close all 

%% Save and view parameters
% set values according to which pieces of data to save or view
view_images = true;
save_localMax_fig = false;
save_ISAR_fig = false;
save_tables = false;
target_folder = 'DSIC_correct/Dataset1';                  % Folder name to save images to

%% Set up parameters
dataset_number = 1;
window_length = 32;                                       % Initial CPTWL                                                              % n value for Martorella's algorithm
hop_size = 1;                                             % hop size for optimal centre profile selection
k_spec = 0.5;                                             % constant for initial threshold
k_DS = 0.1;                                               % constant for Doppler spread acceptable margins

%% Load radar data and define parameters
if dataset_number == 1
    data_file = 'DAP_2010-10-06_18-00-05_002_Umoya_P872_55212';
elseif dataset_number == 2
    data_file = 'DAP_2010-10-06_18-01-54_002_Umoya_P873_50446';
elseif dataset_number == 3
    data_file = 'DAP_2010-10-06_18-14-21_002_Umoya_P872_55336';
elseif dataset_number == 4
    data_file = 'DAP_2010-10-09_06-50-40_008_Umoya_P873_02570';
elseif dataset_number == 5
    data_file = 'DAP_2010-10-09_06-55-26_010_Umoya_P874_03468';
elseif dataset_number == 6
    data_file = 'DAP_2010-10-09_06-58-05_012_Umoya_P874_03085';
else
    error('Invalid dataset number provided')
end

radar_data = radarData(data_file);                          % load data into structure

%% Code
[DS_estimate, centre_profile_idx] = radar_data.DS_estimate_all_frames(window_length, hop_size, k_spec); %Doppler spread at all centre profile indices

% select local maxima image contrast values and corresponding centre
% profiles
min_height = mean(DS_estimate);

% identify local maxima
[DS_local_max, TF] = findpeaks(DS_estimate, "MinPeakHeight", min_height);

% plot local maxima
plot(centre_profile_idx, DS_estimate)
local_max_centre_profile_idx = centre_profile_idx(TF);
hold on
DS_fig = scatter(local_max_centre_profile_idx, DS_local_max);
hold off
xlabel("Centre Profile Index (1)")
ylabel("Doppler Spread (Hz)")
xlim([centre_profile_idx(1) centre_profile_idx(end)])
legend("DS values", "Local Maxima")

% save local maxima figure
if save_localMax_fig
    saveas(DS_fig, [target_folder,'\DSFunction_D',int2str(dataset_number),'.fig']);
end

% Calcualte image contrast at each local DS maxima index
IC_initial = zeros(1, 1);
for i = 1:length(local_max_centre_profile_idx)
    frame = radar_data.select_frame(local_max_centre_profile_idx(i), window_length);
    image = radar_data.ISAR_image_limitedDynamicRange(frame);
    IC_initial(i) = radar_data.image_contrast(image);
end

figure
% select local maxima image contrast values and corresponding centre
% profiles
min_height = mean(IC_initial);

% identify local maxima
[IC_local_max, TF] = findpeaks(IC_initial, "MinPeakHeight", min_height);

% plot local maxima
IC_local_max_centre_profile_idx = local_max_centre_profile_idx(TF);
plot(local_max_centre_profile_idx, IC_initial)
hold on
IC_fig = scatter(IC_local_max_centre_profile_idx, IC_local_max);
hold off
xlabel("Centre Profile Index (1)")
ylabel("Image Contrast (1)")
xlim([local_max_centre_profile_idx(1) local_max_centre_profile_idx(end)])
legend("IC values", "Local Maxima")

if save_localMax_fig
    saveas(IC_fig, [target_folder,'\ICFunction_D',int2str(dataset_number),'.fig']);
end

% Calculate optimal CPTWL at each optimal centre profile index
[DS, IC, window_lengths] = radar_data.DS_IC_window_selection(window_length, IC_local_max_centre_profile_idx, k_DS);

% Identify optimal CPIs with ICs above the 50th percentile of IC values
TF = IC > quantile(IC, 0.5);
optimal_centre_profile_idx = IC_local_max_centre_profile_idx(TF);
optimal_window_lengths = window_lengths(TF);
optimal_IC = IC(TF);

%% Display ISAR images
figure
for i = 1:length(optimal_centre_profile_idx)
    centre_prof = optimal_centre_profile_idx(i);
    window_length = optimal_window_lengths(i);
    frame = radar_data.select_frame(centre_prof, window_length);
    image_limitedDynamicRange_dB = radar_data.ISAR_image_limitedDynamicRange(frame);
    fig = radar_data.plot_ISAR_image(image_limitedDynamicRange_dB);
    title("ISAR image from Dataset "+ dataset_number)
    subtitle("centre profile = "+ centre_prof + ", window length = " + window_length + ", IC = " + optimal_IC(i))
    hold on
    if view_images
        pause(1)
    end
    if save_ISAR_fig
        fileName = [target_folder,'\ISARimage_D',int2str(dataset_number),'_prof',int2str(centre_prof),'_length',int2str(window_length), '_IC', int2str(optimal_IC(i)),'.fig'];
        saveas(fig, fileName)
    end
end

centre_profiles = transpose(optimal_centre_profile_idx);
window_length = transpose(optimal_window_lengths);
IC = transpose(optimal_IC);
T = table(centre_profiles, window_length, IC)
if save_tables
    fileLoc = [target_folder,'\resultsTable_D',int2str(dataset_number)];
    writetable(T, fileLoc);
end

number_of_images = length(optimal_centre_profile_idx);
smallest_window_length = min(window_lengths);
largest_window_length = max(window_lengths);
average_window_length = mean(window_lengths);
smallest_IC = min(IC);
largest_IC = max(IC);
average_IC = mean(IC);
T_sum = table(number_of_images, smallest_window_length, largest_window_length, average_window_length, smallest_IC, largest_IC, average_IC)
if save_tables
    fileLoc = [target_folder,'\summaryTable_D',int2str(dataset_number)];
    writetable(T_sum, fileLoc);
end

beep