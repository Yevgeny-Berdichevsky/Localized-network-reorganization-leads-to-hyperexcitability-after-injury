function [firing_rate] = response_function_arr(input, threshold, slope, max_rate)
%returns firing rate based on piecewise linear function of threshold
%operates on arrays of neurons, returns an array of neurons

firing_rate = max_rate*slope *(input-threshold);
firing_rate(input < threshold) = 0;
firing_rate(input > (threshold + 1/slope)) = max_rate;

end