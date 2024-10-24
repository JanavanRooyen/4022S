classdef radarData
    %radarData Allows for the implementation of CPI selection algorithms on
    %the radar datasets provided by the CSIR.

    properties
        HRR_profiles_all {mustBeNumeric}
        dynamic_range {mustBeNumeric}
        num_of_profiles {mustBeNumeric}
        num_of_range_bins {mustBeNumeric}
        profile_repetition_frequency {mustBeNumeric}
        range_axis {mustBeNumeric}
        frequency_axis {mustBeNumeric}
    end

    methods
        function obj = radarData(data_file)
            %obj Construct an instance of this class
            %   Create object from the input object name
            structure = load(data_file);
            obj.HRR_profiles_all = structure.sb_HRR.G1.HRR_NoMC_calib.';
            obj.dynamic_range = 45;
            obj.num_of_profiles = size(obj.HRR_profiles_all, 1);
            obj.num_of_range_bins = size(obj.HRR_profiles_all, 2);
            obj.profile_repetition_frequency = 1/structure.sb_HRR.G1.Pattern_time;
            obj.range_axis = structure.sb_HRR.G1.xaxis_downrange_m;
            obj.frequency_axis = (-obj.num_of_profiles/2:1:obj.num_of_profiles-1)*obj.profile_repetition_frequency/obj.num_of_profiles;
        end

        function frame = select_frame(obj, centre_profile, window_width)
            %select_frame Output the HRRP data of the specified frame
            %   Calculate the start and end profiles of the specified
            %   window and output the matrix data.
            start_profile = centre_profile - window_width/2;
            end_profile = centre_profile + window_width/2-1;
            frame = obj.HRR_profiles_all(start_profile:end_profile, :);
        end

        function ISAR_image = ISAR_image_complex(obj, frame)
            %ISAR_image_fullDynamicRange Outputs a motion-compensated ISAR
            %image with complex values, without normalised and limited dynamic
            %range
            num_of_profiles = size(frame, 1);
            num_of_range_bins = size(frame, 2);
            frame_aligned = obj.haywood_range_align(frame);                                         % range alignment
            frame_focused = obj.Yuan_phase_adjustment(frame_aligned);                               % phase adjustment
            window_matrix = repmat(hamming(num_of_profiles), 1, num_of_range_bins);
            frame_windowed = frame_focused.*window_matrix;                                          % perform windowing
            ISAR_image = fftshift(fft(frame_windowed, [], 1), 1);                                   % FFT
        end

        function ISAR_image = ISAR_image_limitedDynamicRange(obj, frame)
            %ISAR_image_fullDynamicRange Outputs a motion-compensated ISAR
            %image with dB values, without normalised and limited dynamic
            %range
            num_of_profiles = size(frame, 1);
            num_of_range_bins = size(frame, 2);
            frame_aligned = obj.haywood_range_align(frame);                                         % range alignment
            frame_focused = obj.Yuan_phase_adjustment(frame_aligned);                               % phase adjustment
            window_matrix = repmat(hamming(num_of_profiles), 1, num_of_range_bins);
            frame_windowed = frame_focused.*window_matrix;                                          % perform windowing
            complex_ISAR_image = fftshift(fft(frame_windowed, [], 1), 1);                           % FFT
            ISAR_image = obj.Normalise_limitDynamicRange_ISAR_dB(complex_ISAR_image); % convert image to dB and perform dynamic limiting
        end

        function fig = plot_ISAR_image(obj, ISAR_image_limitedDynamicRange_dB)
            %plot_ISAR_image Plot the heatmap image of the input dB matrix
            %Plot the range-frequency heatmap image of the input dB matrix, 
            %post motion compensation, normalisation and dynamic range 
            %limiting. Output the figure object.
            frame_num_of_profiles = size(ISAR_image_limitedDynamicRange_dB, 1);
            fig = imagesc(obj.range_axis, [-frame_num_of_profiles/2:1:frame_num_of_profiles/2-1]*obj.profile_repetition_frequency/frame_num_of_profiles, ISAR_image_limitedDynamicRange_dB);
            colormap('jet');
            axis xy;
            cb = colorbar;
            title(cb, 'dB')
            xlabel("Range (m)")
            ylabel("Frequency (Hz)")
        end

        function fig = plot_frame_loc(obj, centre_profile_idx, window_width)
            %plot_frame_loc Plots the ISAR image of the provided frame.
            %The specified frame is isolated, converted to an aligned ISAR
            %image and plotted.
            frame = obj.select_frame(centre_profile_idx, window_width);                         % isolate CPI
            image = obj.ISAR_image_limitedDynamicRange(frame);                                  % generate dB image
            fig = obj.plot_ISAR_image(image);                                                   % plot dB image as heatmap
        end

        function [IC, centre_profile_idx] = IC_all_frames(obj, window_length, hop_size)
            %IC_all_frames Outputs an array of image contrast values for
            %all frames in the dataset 
            %Frames of specified length and at the specified spacing are
            %isolated. Motion compensated, normalised ISAR images with
            %limited dynamic ranges are generated, the image contrast of
            %which is calculated. The image contrast, along with the
            %respective centre profile value is outputted.
            frames = floor((obj.num_of_profiles - window_length)/(hop_size) + 1);               % calculate number of frames

            IC = zeros(1, frames);
            centre_profile_idx = zeros(1, frames);

            % iterate through frames, generate ISAR images and compute IC
            for i = 1:frames
                frame = obj.HRR_profiles_all((i-1)*hop_size+1:(i-1)*hop_size+window_length, :); % isolate frame from dataset
                image = obj.ISAR_image_limitedDynamicRange(frame);                              % generate ISAR image of frame
                IC(i) = obj.image_contrast(image);                                                  % Calculate image contrast
                centre_profile_idx(i) = (i-1)*hop_size+1 + window_length/2;                     % determine frame i's centre profile index
            end
        end

        function [optimal_IC, optimal_window_lengths] = IC_window_selection(obj, initial_window_length, centre_profile_idx, IC, n)
            %IC_window_selection Adjust window lengths of specified frame
            %locations to optimise image contrast.
            %Iteratively increase or decrease the window length at each
            %provided centre profile value until the image contrast of the
            %resultant ISAR images stops increasing. Algorithm described [in
            %Martorella, M., & Berizzi, F. (2005). Time Windowing for
            %Highly Focused ISAR Image Reconstruction. IEEE Transations on
            %Aerospace and Electronic Systems, 41(3), 922-1007]
            
            optimal_window_lengths = zeros(1, 1);
            optimal_IC = zeros(1, 1);
            
            % Calculate optimal CPTWL for each centre profile index
            for i = 1:length(centre_profile_idx)
                current_window_length = initial_window_length;
                current_IC = IC(i);
                polarity = 1;
                k = 0;
                first_iteration = true;

                while k < n                                                                         % iterate for k = 0:n-1
                    next_window_length = current_window_length + polarity*2*(n-k);                  % next potential CPTWL
                    
                    % Determine if CPTWL at current centre profile index
                    % exceeds dataset bounds
                    if next_window_length/2 > centre_profile_idx(i)
                        break
                    elseif centre_profile_idx(i) + next_window_length/2-1 > obj.num_of_profiles
                        break
                    elseif next_window_length <= 0
                        break
                    end
                    
                    frame = obj.select_frame(centre_profile_idx(i), next_window_length);            % isolate CPI
                    image = obj.ISAR_image_limitedDynamicRange(frame);                              % generate dB image
                    next_IC = obj.image_contrast(image);                                                % calculate image contrast
                    
                    % set current IC and CPTWL as optimal if IC increased
                    if next_IC > current_IC
                        current_IC = next_IC;
                        current_window_length = next_window_length;
                    
                    elseif next_IC < current_IC && first_iteration
                        polarity = -1;                                                              % decrement for future iterations if IC decreased upon very first iteration
                    
                    else
                        k = k + 1;                                                                  % increment k
                    end
    
                    first_iteration = false;                                                        % false if no longer the first iteration

                end

                % output optimal CPTWL and image contrast for each input
                % centre profile index
                optimal_window_lengths(i) = current_window_length;
                optimal_IC(i) = current_IC;
            end
        end

        function [DS] = doppler_spread(obj, ISAR_image_limitedDynamicRange_dB)
            %doppler_spread calculates the input image's doppler spread
            %according to the optimal threshold.
            %A Doppler spread curve is calculated using potential threshold
            %values. A hyperbola is fitted to the normalised curve, and the
            %vertex calculated. The vertex, once recalibrated to the
            %initial dataset sizes, is considered the optimal threshold
            %according to which the output Doppler spread is calculated.

            num_of_profiles_frame = size(ISAR_image_limitedDynamicRange_dB, 1);
            doppler_spacing = obj.profile_repetition_frequency/num_of_profiles_frame;
            
            profile_sum = sum(ISAR_image_limitedDynamicRange_dB, 2);                           % sum profiles
            thresholds = linspace(min(profile_sum), max(profile_sum), 50);                     % define thresholds array
            DS_all = zeros(1, length(thresholds));

            % calculate initial Doppler spread at each threshold value
            for i = 1:length(thresholds)
                num_profile_sums = length(profile_sum(profile_sum > thresholds(i)));
                DS_all(i) = num_profile_sums*doppler_spacing;
            end
            
            %Normalise between 0 and 1
            x_data = (thresholds - min(profile_sum))/abs(max(profile_sum) - min(profile_sum));
            y_data = DS_all/DS_all(1);

            %% fit hyperbola to data
            initial_guess = [7, 0.8];                                                           % initial values for parameters a and b in hyperbola function
            hyperbola = @(params, x) 1 ./ (params(1)*x+params(2));                              % hyperbola function
            opt = optimoptions('lsqcurvefit', 'Display', 'off', 'Algorithm', 'levenberg-marquardt');
            [params_opt] = lsqcurvefit(hyperbola, initial_guess, x_data, y_data, [], [], [], [], [], [], [], opt);                 % optimise parameters

            % Extract optimized parameters
            a = params_opt(1);
            b = params_opt(2);

            %% Calculate optimal threshold value
            vertex = (sqrt(a) - b)/a;                                                           % hyperbola's vertex value
            threshold_opt = min(profile_sum) + vertex * (max(profile_sum) - min(profile_sum));  % reverse normalisation of vertex value
            
            %% Calculate final doppler spread according to 
            DS = length(profile_sum(profile_sum > threshold_opt)) * doppler_spacing;
        end

        function [DS, centre_profile_idx] = DS_estimate_all_frames(obj, window_length, hop_size, k_spec)
            %DS_estimate_all_frames estimates the Doppler spread, with a
            %non-optimal threshold, for all CPIs within the dataset.
            % A spectrogram of the entire dataset with summed profiles is
            % used to estimate the Doppler spread.
            
            %apply motion compensation algorithms
            HRRP_aligned = obj.haywood_range_align(obj.HRR_profiles_all);
            HRRP_focused = obj.Yuan_phase_adjustment(HRRP_aligned);

            %% Spectrogram
            % define parameters
            frames = floor((obj.num_of_profiles - window_length)/(hop_size) + 1);           % calculate number of frames
            doppler_spacing = obj.profile_repetition_frequency/obj.num_of_profiles;         % define frequency resolution

            DS = zeros(1, 1);                                                               % array to store DS estimates of each frame's image
            centre_profile_idx = zeros(1, 1);                                               % array to store centre profile indexes
            w = hann(window_length);                                                        % window to apply in spectrogram calculation

            HRRP_sum = sum(HRRP_focused, 2);

            % compute DS of each column of the spectrogram
            for i = 1:frames
                % apply windowing
                frame_sum = [(HRRP_sum(1:(i-1)*hop_size))*0; HRRP_sum((i-1)*hop_size+1:(i-1)*hop_size+window_length).*w; HRRP_sum((i-1)*hop_size+window_length+1:obj.num_of_profiles)*0];
                fft_sum = 20*log10(abs(fft(frame_sum, [], 1)));                             % calculate FFT of summed profiles
                centre_profile_idx(i) = (i-1)*hop_size+1 + window_length/2;                 % determine frame i's centre profile index
                
                % statistical values of fft_profile_sum
                sigma = std(fft_sum);
                mu = mean(fft_sum);
                lambda = mu-k_spec*sigma;
                
                DS(i) = sum(fft_sum<lambda)*doppler_spacing;                                % determine range of doppler frequency values from number of elements within the specified range
            end
        end

        function [optimal_DS, optimal_IC, optimal_window_lengths] = DS_IC_window_selection(obj, initial_window_length, centre_profile_idx, k_DS)        
        %DS_IC_window_selection Calcualates the optimal CPTWL for all the
        %provided centre profile indices.
        %For each input centre profile index, the control Doppler spread
        %and appropriate margin for error calculated. Thereafter, a range
        %of potential CPTWLs is tested, and the value that produces an
        %image with a Doppler spread within defined tolerances and
        %producing the maximum image contrast is selected as the optimal
        %CPTWL for that centre profile index. An array of optimal CPTWLs,
        %corresponding to the input profiles, is outputted, along with the
        %corresponding Doppler spread and image contrast.

        optimal_DS = zeros(1,length(centre_profile_idx));
        optimal_window_lengths = zeros(1,length(centre_profile_idx));
        optimal_IC = zeros(1,length(centre_profile_idx));

        % perform CPTWL selection for each input profile index
        for i = 1:length(centre_profile_idx)
            frame = obj.select_frame(centre_profile_idx(i), initial_window_length);             %isolate frame
            image = obj.ISAR_image_limitedDynamicRange(frame);                                  % generate dB image
            DS_control = obj.doppler_spread(image);                                             % control, accurate Doppler spread
            tolerance_lower = DS_control - k_DS*DS_control;                                     % define tolerances
            tolerance_upper = DS_control + k_DS*DS_control;

            window_length = initial_window_length;
            IC = 0;
            
            lower_bound = centre_profile_idx(i)*2-1;
            upper_bound = 2*(obj.num_of_profiles + 1 - centre_profile_idx(i));
            potential_window_lengths = min(min(32, lower_bound), upper_bound):2:min(min(lower_bound, upper_bound), 200);    % CPTWLs to be tested
            for j = 1:length(potential_window_lengths)                                          % test each potential CPTWL
                next_window_length = potential_window_lengths(j);
                
                frame = obj.select_frame(centre_profile_idx(i), next_window_length);            % new CPI with potential CPTWL
                image = obj.ISAR_image_limitedDynamicRange(frame);                              % generate dB image
                next_DS = obj.doppler_spread(image);                                            % calculate new accurate Doppler spread

                if next_DS < tolerance_upper && next_DS > tolerance_lower                       % test if value falls within tolerances
                    next_IC = obj.image_contrast(image);                                            % calculate image contrast of new image
                    if next_IC > IC                                                             % if image contrast increased, store current
                        IC = next_IC;                                                           % image contrast, Doppler spread and CPTWL as optimal
                        window_length = next_window_length;
                        DS = next_DS;
                    end
                end
            end

            % optimal output values for each centre profile index
            optimal_DS(i) = DS;
            optimal_IC(i) = IC;
            optimal_window_lengths(i) = window_length;

        end

        end

        function [range_aligned_image] = haywood_range_align(obj, frame_complex)
        %range_alignment: aligns the range profiles of the input ISAR image
        %   The profile with the dominant scatterer closest to the image centre is
        %   defined as the reference profile. All profiles in the input image are 
        %   cross-correlated with the reference profile to calculate their required
        %   delay. Said delay is then fitted to a low-order polynomial and then
        %   applied to the input image through phase adjustment in the frequency
        %   domain. The output is the range-aligned ISAR image.             
        
        % define CPI parameters
        num_of_profiles = size(frame_complex, 1);
        num_of_range_bins = size(frame_complex, 2);

        centre = (num_of_range_bins+1)/2;                                      % define centre range bin
        
        %find profile with dominant scatterer nearest to centre range bin
        [~, peak_range] = max(frame_complex, [], 2);                            % find range bin of dominant scatterers in each profile
        bins_from_centre = abs(centre - peak_range);                            % distance between centre range bin and dominant scatterer range bins
        [~, closest_profile] = min(bins_from_centre);                           % number of profile with smallest distance between dominant scatterer and centre range bins
        ref_profile = frame_complex(closest_profile,:);                         % define reference profile with obtained profile index
        
        %perform cross correlation to obtain delay
        corr = ifft(fft(abs(frame_complex), [ ], 2).*repmat(conj(fft(abs(ref_profile))), num_of_profiles, 1), [], 2);
        
        %calculate delays for each profile  
        [~, delay] = max(corr, [], 2);                                          % find range bins of dominant scatterers, and therefore required delay
        delay = delay - 1;                                                      % adjust delay values to fall within [0; n-1]
        
        %define phi matrix
        m = 0:1:(num_of_range_bins-1);
        delay_unwrapped = unwrap(delay, (num_of_range_bins-1));                % unwrap delay function to eliminate jumps
        [p, ~, mu] = polyfit(1:num_of_profiles, delay_unwrapped, 1);           % fit polynomial to the unwrapped delays
        d = transpose(polyval(p, 1:num_of_profiles, [], mu));                  % evaluate the value of the polynomial for each profile
        phi = [exp((1j*2*pi/num_of_range_bins)*(d*m))];                        % calculate phi matrix, different polynomial value in each row and m value in each column
    
        range_aligned_image = ifft(phi.*fft(frame_complex, [], 2), [], 2);      % shift image with phi matrix

        end

        function [phase_adjusted_image] = Yuan_phase_adjustment(obj, range_aligned_image_complex)
        %Yuan_phase_adjustment performs weighted multiple scatterer phase
        %adjustment to the input range-aligned Range-profile image.
        %   Stable scatterers are chosen according to the criteria 
        %   var/(var + mean^2), which selects for prominent scatterers that vary 
        %   minimally between pulses. The first 11 scatterers with the highest
        %   stability (smallest calculated values) are chosen. For each complex 
        %   scatterer value in each profile, it is multipled by the conjugate of
        %   the same scatterer's value in the first profile. The mean complex value
        %   of all products in each profile is calculated, the phase of which is
        %   retained. This phase is used to form a complex exponential correction
        %   vector to apply to the input image.
        
        %% Identify scatterers
        range_bin_mean = mean(abs(range_aligned_image_complex), 1);                                                 % calculate mean magnitude of each range bin
        range_bin_var = var(abs(range_aligned_image_complex), [], 1);                                               % calculate magnitude variance of each range bin
        stability = (range_bin_var)./((range_bin_var) + (range_bin_mean.^2));                                       % critera for low variance and high mean (i.e. high stability)
        
        [stability_ascending, range_bin] = sort(stability);                                                         % sort stability values in ascending order
        scatterer_range_bins = range_bin(1:min([find(stability_ascending < 0.16, 1, 'last'),  11]));                % select all range bins, capped at first 11, with a var_mean_ratio < 0.16
        
        %% Calculate and apply phase correction
        scatterers = range_aligned_image_complex(:,scatterer_range_bins);                                           % all complex values of scatterers
        
        profile_zero_mat = repmat(conj(scatterers(1,:)), size(range_aligned_image_complex, 1), 1);                  % create matrix of profile 0's conjugate
        phi = angle(mean((scatterers).*profile_zero_mat, 2));                                                       % phase of mean product of conjugate of profile 0 and profile N-1
        phase_correction_vector = exp(-1j*phi);                                                                     % complex exponential with the phase of -phi
        phase_adjusted_image = range_aligned_image_complex.*phase_correction_vector;                                % apply phase adjustment vector to each respective profile of the image
        
        end
        
        function [Normalise_limitDynamicRange_ISAR_dB] = Normalise_limitDynamicRange_ISAR_dB(obj, ISAR_image_complex_linear)
        %Normalise_limitDynamicRange_ISAR_dB converts the input, complex
        %valued ISAR image to dB and performs dynamic range limiting.
        
        ISAR_image_complex_linear = abs(ISAR_image_complex_linear)./max(max(abs(ISAR_image_complex_linear)));       % normalise values
        ISAR_image_dB = 20*log10(abs(ISAR_image_complex_linear));                                                   % convert to dB
        max_value_inplot = max(max(ISAR_image_dB));                                                                 % find maximum dB value in image
        indx = find(ISAR_image_dB < (max_value_inplot-obj.dynamic_range));                                          % find index of values smaller than the maximum value minus the dynamic range limitation
        ISAR_image_dB(indx) = (max_value_inplot-obj.dynamic_range);                                                 % set value at indices to difference between the maximum value and the dynamic range limitation
        ISAR_image_dB =  ISAR_image_dB - max_value_inplot;                                                          % subtract maximum value from all points in the image
          
        Normalise_limitDynamicRange_ISAR_dB = ISAR_image_dB;
        end

        function [IC] = image_contrast(obj, ISAR_image_dB)
        %image_contrast calculates the IC value for a given ISAR image that's
        %already been windowed appropriately.
        %   The intensity of each datum in the ISAR image is obtained. The image 
        %   contrast is calculated using both each intensity value as well as the 
        %   average of all the intensity values of the image.
        
        intensity = (10.^(ISAR_image_dB./20)).^2;                                                               % intensity of each datum of entire image
        
        mean_intensity = mean(intensity, "all");                                                                % average of entire image's intensity
        IC = sqrt(mean((intensity - mean_intensity).^2, "all"))/mean_intensity;                                 % image contrast

        end

    end
end