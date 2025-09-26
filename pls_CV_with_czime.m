
function mdl = pls_CV(y,x,npls_comps, n_outer_reps,n_inner_reps)
y = els_high_matrix
x = [ADRA1A, ADRA1B, ADRA1D,ADRA2A,ADRA2B,ADRA2C,ADRB1,ADRB2,ADRB3,CHRM1,CHRM2,CHRM3,CHRM4,CHRM5,CHRNA2,CHRNA3,CHRNA4,CHRNA5,CHRNA6,CHRNA7,CHRNA9,CHRNA10,CHRNB2,CHRNB3,D1,D2,D3,D4,D5,HT1A,HT1B,HT1B,HT1D,HT1E,HT1F,HT2A,HT2B,HT2C,HT3A,HT4,HT5A,HT6,HT7]
npls_comps = 3
n_outer_reps = 50
n_inner_reps = 2

% function for pls regression with nested k-fold CV employed to tune number
% of components and validate the model (note usage of coefficient of
% determination, other metrics such as correlation and mean squared error
% are also suitable)
%
% script will: define CV folds, divide data into train/test groups, start
% inner loop to further divide training data into inner_train/inner_test, use that to
% find out how many components are optimal, once best number of components
% is found, will use that number (minimum of 1) to fit the model to the
% whole data and make predictions of outter test sample
%
% y - matrix of subject*outcomes that will be predicted by plsr
% x - matrix of subject*predictors that will be used to predict y
% npls_comps - max number of pls components explored in inner cv - note: it
% should be redundant to explore large numbers of components (i.e.
% approaching to npls_comps = MIN(SIZE(X,1)-1, SIZE(X,2)))
% n_outer_reps - number of outer CV repeats for model validation
% n_inner_reps - number of outer CV repeats for hyperparameter tuning

nSubs = length(y);
% nested crossvalidation for hyperparameter (NCOMPS) tuning
k_outer = 10;
k_inner = 5;

for i_outer = 1:n_outer_reps

    % define cv folds
    cv_indices_outer = crossvalind('Kfold', nSubs, k_outer); % on each repetition creates a new index

    for j_outer = 1:k_outer % for each cross validation fold

        % divide training and validation folds
        test = (cv_indices_outer == j_outer);
        train = ~test;

        x_train = x(train,:);
        y_train = y(train,:);

        x_test = x(test,:);
        y_test = y(test,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% tuning - optimal number of pls components %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        for i_inner = 1:n_inner_reps
            % define inner cv folds
            cv_indices_inner = crossvalind('Kfold', length(y_train), k_inner);
            for j_inner = 1:k_inner % for each cross validation fold

                % divide training and validation folds
                test_inner = (cv_indices_inner == j_inner);
                train_inner = ~test_inner;

                x_train_inner = x_train(train_inner,:);
                y_train_inner = y_train(train_inner,:);

                x_test_inner = x_train(test_inner,:);
                y_test_inner = y_train(test_inner,:);

                for comps = 1:npls_comps % predict unseen scores with each n_pls_comps

                    % fit pls regression with target
                    % normalize here x_train_inner and y_train_inner

                    % normalize data
                    [Z, mean_x_train_inner, standard_deviation_x_train_inner] = zscore(x_train_inner);
                    % in the presence of 0-connections obtained due to thresholding
                    Z(isinf(Z)) =[];
                    Z(isnan(Z)) =[];
                    % repeat for response variable (y)
                    [z, mean_y_train_inner, standard_deviation_y_train_inner] = zscore(y_train_inner);
                    z(isinf(z)) =[];
                    z(isnan(z)) =[];

                    [~,~,~,~,BETA,~,~] = plsregress(Z, z, comps);

                    % apply the normalisation parameters (mean and standard
                    % deviation) to the  x_test_inner and y_test_inner

                    x_test_inner = bsxfun(@minus,x_test,mean_x_train_inner);
                    x_test_inner = bsxfun(@rdivide,x_test,standard_deviation_x_train_inner);

                    y_test_inner = bsxfun(@minus,x_test,mean_y_train_inner);
                    y_test_inner = bsxfun(@rdivide,x_test,standard_deviation_y_train_inner);



                    % predict unseen sample data, BETA includes constant
                    ests_inner = [ones(size(x_test_inner,1),1) x_test_inner]*BETA;

                    SStot = sum((y_test_inner - mean(y_test_inner)).^2);
                    SSres = sum((y_test_inner - ests_inner).^2);
                    e = 1 - SSres./SStot; % coefficient of determination

                    % e = mean((y_test_inner - ests_inner).^2);% Mean Squared Error
                    % e = corr(y,x);

                    err(j_inner,comps) = mean(e);
                    clear Z z BETA
                end % end inner - number of comps to consider
                err_cv= mean(err);
            end % end inner - CV folds
        end % end inner - CV repeats
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% fit the model using n-comps that are needed to produce minimum
        %%%%%%%% error during validation of tuning sample %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [M,I] = max(mean(err_cv)); % find number of components with minimum error

        % normalize here x_train and y_train

        % normalize data
        [Z, mean_x_train, standard_deviation_x_train] = zscore(x_train);
        % in the presence of 0-connections obtained due to thresholding
        Z(isinf(Z)) =[];
        Z(isnan(Z)) =[];
        % repeat for response variable (y)
        [z, mean_y_train, standard_deviation_y_train] = zscore(y_train);
        z(isinf(z)) =[];
        z(isnan(z)) =[];
        % normalize here x_test and y_test

        x_test = bsxfun(@minus,x_test,mean_x_train);
        x_test = bsxfun(@rdivide,x_test,standard_deviation_x_train);

        y_test = bsxfun(@minus,x_test,mean_y_train);
        y_test = bsxfun(@rdivide,x_test,standard_deviation_y_train);

        % fit pls model using that many components
        [~,~,~,~,BETA,PCTVAR,~,STATS] = plsregress(Z, z, I);

        % predict scores of validation fold
        ests = [ones(size(x_test,1),1) x_test]*BETA;
        err_outer = mean((y_test - ests).^2);


        SStot = sum((y_test - mean(y_test)).^2);
        SSres = sum((y_test - ests).^2);
        CoD = 1 - SSres./SStot;



        % save what is needed
        mdl.n_comps(i_outer,j_outer) = I; % number of components
        mdl.css_min(i_outer,j_outer) = M; % value
        %out.css.L{i,j} = comp_L(:,I); % loss at this value
        mdl.pls.BETA{i_outer,j_outer} = BETA; % BETA is a (p + 1)-by-m matrix, where p
        %is the number of x data columns prior to feature selection
        mdl.mse_outer{i_outer,j_outer} = err_outer;
        mdl.R_outer{i_outer,j_outer} = CoD;
        
    end
end
end
