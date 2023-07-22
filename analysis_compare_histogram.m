%histogram 1
load('excursion_result_time_poisson_different_lambda22-Apr-2023 21_45_34.mat')

%mean_time_to_excursion for 10000 trails
list_of_ones = find(excursion_result == 1);
mean_time_to_excursion = sum(excursion_time(list_of_ones,:)./length(list_of_ones))

%histogram
subplot(2,3,1)
histogram_6 = histogram(excursion_time(list_of_ones), 30)
xlabel('mean time to excursion')
ylabel('number of trails')
title('mean of \tau = 1')
%lambda = 600



%histogram 2
load('excursion_result_time_poisson_lambda14-Apr-2023 14_18_56.mat')
%this loads three variables: excursion_time, excursion_result, which_lambda for
%5000 trials
%before we load the next file, we rename the variables so they don't get
%overwritten
excursion_time1 = excursion_time;
excursion_result1 = excursion_result;
which_lambda1 = which_lambda;

%load the next file 5000 trails
load('excursion_result_time_poisson_different_lambda18-Apr-2023 13_11_51.mat')
excursion_time2 = excursion_time;
excursion_result2 = excursion_result;
which_lambda2 = which_lambda;


%combine the results of variables in the three files
excursion_result = [excursion_result1; excursion_result2];
excursion_time = [excursion_time1; excursion_time2];
which_lambda = [which_lambda1; which_lambda2];
%lambda = 85.7143

%mean_time_to_excursion for 10000 trails
list_of_ones = find(excursion_result == 1);
mean_time_to_excursion = sum(excursion_time(list_of_ones,:)./length(list_of_ones))

%histogram
%figure
subplot(2,3,2)
histogram_1 = histogram(excursion_time(list_of_ones), 30)
xlabel('mean time to excursion')
ylabel('number of trails')
title('mean of \tau = 7')



%histogram 3
load('excursion_result_time_poisson_different_lambda18-Apr-2023 14_15_49.mat')
%this loads three variables: excursion_time, excursion_result, which_lambda for
%5000 trials
%before we load the next file, we rename the variables so they don't get
%overwritten
excursion_time1 = excursion_time;
excursion_result1 = excursion_result;
which_lambda1 = which_lambda;

%load the next file 5000 trails
load('excursion_result_time_poisson_lambda14-Apr-2023 15_23_19.mat')
excursion_time2 = excursion_time;
excursion_result2 = excursion_result;
which_lambda2 = which_lambda;

%combine the results of variables in the three files
excursion_result = [excursion_result1; excursion_result2];
excursion_time = [excursion_time1; excursion_time2];
which_lambda = [which_lambda1; which_lambda2];
%lambda = 60

%mean_time_to_excursion for 10000 trails
list_of_ones = find(excursion_result == 1);
mean_time_to_excursion = sum(excursion_time(list_of_ones,:)./length(list_of_ones))

%histogram
subplot(2,3,3)
histogram_2 = histogram(excursion_time(list_of_ones), 30)
xlabel('mean time to excursion')
ylabel('number of trails')
title('mean of \tau = 10')



%histogram 4
load('excursion_result_time_poisson_different_lambda18-Apr-2023 15_12_40.mat')
%this loads three variables: excursion_time, excursion_result, which_lambda for
%5000 trials
%before we load the next file, we rename the variables so they don't get
%overwritten
excursion_time1 = excursion_time;
excursion_result1 = excursion_result;
which_lambda1 = which_lambda;

%load the next file 5000 trails
load('excursion_result_time_poisson_lambda14-Apr-2023 16_23_02.mat')
excursion_time2 = excursion_time;
excursion_result2 = excursion_result;
which_lambda2 = which_lambda;

%combine the results of variables in the three files
excursion_result = [excursion_result1; excursion_result2];
excursion_time = [excursion_time1; excursion_time2];
which_lambda = [which_lambda1; which_lambda2];
%lambda = 46.1538

%mean_time_to_excursion for 10000 trails
list_of_ones = find(excursion_result == 1);
mean_time_to_excursion = sum(excursion_time(list_of_ones,:)./length(list_of_ones))

%histogram
subplot(2,3,4)
histogram_3 = histogram(excursion_time(list_of_ones), 30)
xlabel('mean time to excursion')
ylabel('number of trails')
title('mean of \tau = 13')



%histogram 5
load('excursion_result_time_poisson_different_lambda18-Apr-2023 16_05_29.mat')
%this loads three variables: excursion_time, excursion_result, which_lambda for
%5000 trials
%before we load the next file, we rename the variables so they don't get
%overwritten
excursion_time1 = excursion_time;
excursion_result1 = excursion_result;
which_lambda1 = which_lambda;

%load the next file 5000 trails
load('excursion_result_time_poisson_lambda14-Apr-2023 17_17_12.mat')
excursion_time2 = excursion_time;
excursion_result2 = excursion_result;
which_lambda2 = which_lambda;

%combine the results of variables in the three files
excursion_result = [excursion_result1; excursion_result2];
excursion_time = [excursion_time1; excursion_time2];
which_lambda = [which_lambda1; which_lambda2];
%lambda = 37.50

%mean_time_to_excursion for 10000 trails
list_of_ones = find(excursion_result == 1);
mean_time_to_excursion = sum(excursion_time(list_of_ones,:)./length(list_of_ones))

%histogram
subplot(2,3,5)
histogram_4 = histogram(excursion_time(list_of_ones), 30)
xlabel('mean time to excursion')
ylabel('number of trails')
title('mean of \tau = 16')



%histogram 6
load('excursion_result_time_poisson_different_lambda19-Apr-2023 11_59_52.mat')
%this loads three variables: excursion_time, excursion_result, which_lambda for
%5000 trials
%before we load the next file, we rename the variables so they don't get
%overwritten
excursion_time1 = excursion_time;
excursion_result1 = excursion_result;
which_lambda1 = which_lambda;

%load the next file 5000 trails
load('excursion_result_time_poisson_lambda15-Apr-2023 12_51_33.mat')
excursion_time2 = excursion_time;
excursion_result2 = excursion_result;
which_lambda2 = which_lambda;

%combine the results of variables in the three files
excursion_result = [excursion_result1; excursion_result2];
excursion_time = [excursion_time1; excursion_time2];
which_lambda = [which_lambda1; which_lambda2];
%lambda = 6

%mean_time_to_excursion for 10000 trails
list_of_ones = find(excursion_result == 1);
mean_time_to_excursion = sum(excursion_time(list_of_ones,:)./length(list_of_ones))

%histogram
subplot(2,3,6)
histogram_5 = histogram(excursion_time(list_of_ones), 30)
xlabel('mean time to excursion')
ylabel('number of trails')
title('mean of \tau = 100')

