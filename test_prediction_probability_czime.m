
y_fitted = % the predictions of y given by your model
y_observed = % true, observed, measured y that you collected from the pps

obs_loss = % your coeeficient of determination/correlation from the mdl variable
s = sort(obs_loss);

ci_L = s(250); % confidence intervals
ci_U = s(9750);
CI = [ci_L ci_U]; % confidence intervals

p = test_sig(y_fitted, y_observed);

%% another way to do this, edit at own demise
alpha = 0.001;
alpha_lower = alpha / 2;  
alpha_upper = 1 - alpha / 2;  

% Calculate the critical percentiles of the null distribution
lower_percentile = prctile(cell2mat(null_t_values_2back), alpha_lower * 100, 1);  % 0.05th percentile
upper_percentile = prctile(cell2mat(null_t_values_2back), alpha_upper * 100, 1); % 99.95th percentile

% Initialize matrix to store significant results
significant_t_values_2back = zeros(size(observed_t_values));

% Determine significant t-values
for net1 = 1
    for net2 = 1:length(observed_t_values)
        if observed_t_values(net1, net2) < lower_percentile(:, net2) || observed_t_values(net1, net2) > upper_percentile(:, net2)
            significant_t_values_2back(net1, net2) = observed_t_values(net1, net2);
        end
    end
end


%%
function p = test_sig(y_observed, y_fitted) % I deleted the input variable obs_loss, which will be computed inside now

obs_loss = est_loss(obs,predicted(:)); % Computes the loss on the un-randomised data.

predictedv = predicted(:); % This is just you don't have to do this vectorisation 10000 times inside the loop, it will save time.

for b = 1:10000

    rand_obs = obs(randperm(length(obs)));

    rand_loss(b) = est_loss(rand_obs,predictedv);

end

p = (sum(rand_loss>=obs_loss))/length(rand_loss(:));
 
end

 

function out = est_loss(y_test,ests)    

%%% coef of determination

SStot = sum((y_test - mean(y_test)).^2);

SSres = sum((y_test - ests).^2);

out = 1 - SSres./SStot;

%%% pearson's correlation
% out = corr(y_test, ests);

%%% root mean squared error
% out = sqrt(mean((y_test - ests).^2));

end 
