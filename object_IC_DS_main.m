%% Code to save ISAR images identified with Cataldo's IC and DS algorithm, as well as IC function associated with dataset.

clear all
clc
close all 


%% Fill in details
dataset_number = 1;
target_folder = 'DSIC_correct/Dataset1';                  % Folder name to save images to
fileName_suffix = 'h1_kspec05_kds01';
window_length = 32;                                                         % Initial window length for Martorella's algorithm                                                              % n value for Martorella's algorithm
hop_size = 1;
k_spec = 0.5;
k_DS = 0.1;
view_images = false;
save_localMax_fig = false;
save_ISAR_fig = false;
save_tables = false;

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

radar_data = radarData(data_file);

tStart = tic;
%% Functions
[DS_estimate, centre_profile_idx] = radar_data.DS_estimate_all_frames(window_length, hop_size, k_spec);
[DS_local_max, local_max_centre_profile_idx, DS_fig] = local_maxima_DS(DS_estimate, centre_profile_idx);

legend("DS values", "Local Maxima")
if save_localMax_fig
    saveas(DS_fig, [target_folder,'\DSFunction_D',int2str(dataset_number),'_', fileName_suffix,'.fig']);
end

IC_initial = zeros(1, 1);
for i = 1:length(local_max_centre_profile_idx)
    frame = radar_data.select_frame(local_max_centre_profile_idx(i), window_length);
    image = radar_data.ISAR_image_limitedDynamicRange(frame);
    IC_initial(i) = image_contrast(image);
end

figure
[IC_local_max, IC_local_max_centre_profile_idx] = local_maxima_DS(IC_initial, local_max_centre_profile_idx);

legend("IC values", "Local Maxima")
if save_localMax_fig
    saveas(DS_fig, [target_folder,'\ICFunction_D',int2str(dataset_number),'_', fileName_suffix,'.fig']);
end

[DS, IC, window_lengths] = radar_data.DS_IC_window_selection(window_length, IC_local_max_centre_profile_idx, k_DS);

% Identify collection of optimal frames based on IC
figure
plot(IC_local_max_centre_profile_idx, IC)
%[optimal_IC, TF] = findpeaks(IC, "MinPeakHeight", quantile(IC, 0.8));
[optimal_IC, TF] = findpeaks(IC);
optimal_centre_profile_idx = IC_local_max_centre_profile_idx(TF);
optimal_DS = DS(TF);
optimal_window_lengths = window_lengths(TF);
hold on
IC_fig = scatter(optimal_centre_profile_idx, optimal_IC);
hold off
xlabel("Centre Profile Index (1)")
ylabel("Image Contrast (1)")
xlim([IC_local_max_centre_profile_idx(1) IC_local_max_centre_profile_idx(end)])

ent = zeros(1, 1);

length(IC)
TF = IC > quantile(IC, 0.5);
optimal_centre_profile_idx = IC_local_max_centre_profile_idx(TF);
optimal_window_lengths = window_lengths(TF);
optimal_IC = IC(TF);

t_elapsed = toc(tStart);
%beep
%pause(5)
figure
for i = 1:length(optimal_centre_profile_idx)
    centre_prof = optimal_centre_profile_idx(i);
    window_length = optimal_window_lengths(i);
    frame = radar_data.select_frame(centre_prof, window_length);
    image_limitedDynamicRange_dB = radar_data.ISAR_image_limitedDynamicRange(frame);
    ent(i) = entropy(rescale(image_limitedDynamicRange_dB));
    fig = radar_data.plot_ISAR_image(image_limitedDynamicRange_dB);
    title("ISAR image from Dataset "+ dataset_number)
    subtitle("centre profile = "+ centre_prof + ", window length = " + window_length + ", IC = " + optimal_IC(i) + ", entropy = "+ ent(i))
    hold on
    if view_images
        pause(1)
    end
    if save_ISAR_fig
        fileName = [target_folder,'\ISARimage_D',int2str(dataset_number),'_prof',int2str(centre_prof),'_length',int2str(window_length), '_IC', int2str(optimal_IC(i)), '_', fileName_suffix,'.fig'];
        saveas(fig, fileName)
    end
end

centre_profiles = transpose(optimal_centre_profile_idx);
window_length = transpose(optimal_window_lengths);
IC = transpose(optimal_IC);
image_entropy = transpose(ent);
T = table(centre_profiles, window_length, IC, image_entropy)
if save_tables
    fileLoc = [target_folder,'\resultsTable_Latex_D',int2str(dataset_number),'_',fileName_suffix];
    table2latex(T, fileLoc);
    fileLoc = [target_folder,'\resultsTable_D',int2str(dataset_number),'_',fileName_suffix];
    writetable(T, fileLoc);
end

number_of_images = length(optimal_centre_profile_idx);
smallest_window_length = min(window_lengths);
largest_window_length = max(window_lengths);
average_window_length = mean(window_lengths);
smallest_IC = min(IC);
largest_IC = max(IC);
average_IC = mean(IC);
smallest_ent = min(ent);
largest_ent = max(ent);
average_ent = mean(ent);
T_sum = table(number_of_images, smallest_window_length, largest_window_length, average_window_length, smallest_IC, largest_IC, average_IC, smallest_ent, largest_ent, average_ent, t_elapsed)
if save_tables
    fileLoc = [target_folder,'\summaryTable_Latex_D',int2str(dataset_number),'_',fileName_suffix];
    table2latex(T_sum, fileLoc);
    fileLoc = [target_folder,'\summaryTable_D',int2str(dataset_number),'_',fileName_suffix];
    writetable(T_sum, fileLoc);
end

beep