%% Maximum Contrast based Automatic Time Window Selection
% adjust parameters according to which data to save, and where to.

clear all
clc
close all 

%% Save and view parameters
% set values according to which pieces of data to save or view
view_images = true;
save_IC_fig = false;
save_ISAR_fig = false;
save_tables = false;
target_folder = 'results\Dataset6';                  % Folder name to save images to

%% Set up parameters
dataset_number = 1;
window_length = 64;                                  % initial CPTWL
n = 6;                                               % n value for CPTWL selection
hop_size = 1;                                        % hop size for optimal centre profile selection

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

radar_data = radarData(data_file);                      % load data into structure

%% Code
[IC, centre_profile_idx] = radar_data.IC_all_frames(window_length, hop_size);       % image contrast at all centre profile indices

% select local maxima image contrast values and corresponding centre
% profiles
min_height = quantile(IC, 0.8);
min_dist = 18;

% identify local maxima
[IC_local_max, TF] = findpeaks(IC, "MinPeakHeight", min_height, "MinPeakDistance", min_dist);

% plot local maxima
plot(centre_profile_idx, IC)
optimal_centre_profile_idx = centre_profile_idx(TF);
hold on
IC_fig = scatter(optimal_centre_profile_idx, IC_local_max);
hold off
xlabel("Centre Profile Index (1)")
ylabel("Image Contrast (1)")
xlim([centre_profile_idx(1) centre_profile_idx(end)])
legend("IC values", "Local Maxima")

% save local maxima figure
if save_IC_fig
    saveas(IC_fig, [target_folder,'\ICFunction_D',int2str(dataset_number),'.fig']);
end

% Select optimal CPTWL for each local maxima's centre profile index
[optimal_IC, optimal_window_lengths] = radar_data.IC_window_selection(window_length, optimal_centre_profile_idx, IC_local_max, n);


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

%% Display and save tables
centre_profiles = transpose(optimal_centre_profile_idx);
window_length = transpose(optimal_window_lengths);
IC = transpose(optimal_IC);
T = table(centre_profiles, window_length, IC)
if save_tables
    fileLoc = [target_folder,'\resultsTable_D',int2str(dataset_number)];
    writetable(T, fileLoc);
end

number_of_images = length(optimal_centre_profile_idx);
smallest_window_length = min(optimal_window_lengths);
largest_window_length = max(optimal_window_lengths);
average_window_length = mean(optimal_window_lengths);
smallest_IC = min(optimal_IC);
largest_IC = max(optimal_IC);
average_IC = mean(optimal_IC);
T_sum = table(number_of_images, smallest_window_length, largest_window_length, average_window_length, smallest_IC, largest_IC, average_IC)
if save_tables
    fileLoc = [target_folder,'\summaryTable_D',int2str(dataset_number)];
    writetable(T_sum, fileLoc);
end

beep