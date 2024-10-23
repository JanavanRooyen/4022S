function [IC] = image_contrast(ISAR_image_dB)
%image_contrast calculates the IC value for a given ISAR image that's
%already been windowed appropriately.
%   The intensity of each datum in the ISAR image is obtained. The image 
%   contrast is calculated using both each intensity value as well as the 
%   average of all the intensity values of the image.

intensity = (10.^(ISAR_image_dB./20)).^2;                                                     % intensity of each datum of entire image

mean_intensity = mean(intensity, "all");                                                % average of entire image's intensity
IC = sqrt(mean((intensity - mean_intensity).^2, "all"))/mean_intensity;                 % image contrast

end