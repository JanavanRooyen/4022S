function [local_max,tau_optimal, IC_fig] = local_maxima(func,centre_profile_idx)
%local_maxima finds and local maxima of the provided function. The function
%along with these maxima are plotted.
%   The local maxima values of the provided function are identified
%   according to the parameters set. The function is plotted and its local
%   maxima indicated. These local maxima and their corresponding tau values
%   are returned.

min_height = quantile(func, 0.8);
min_dist = 18;

%% Obtain and plot tau values at local DS maximas
plot(centre_profile_idx, func)
[local_max, TF] = findpeaks(func, "MinPeakHeight", min_height, "MinPeakDistance", min_dist);
tau_optimal = centre_profile_idx(TF);
hold on
IC_fig = scatter(tau_optimal, local_max);
hold off
xlabel("Centre Profile Index (1)")
ylabel("Image Contrast (1)")
xlim([centre_profile_idx(1) centre_profile_idx(end)])

end