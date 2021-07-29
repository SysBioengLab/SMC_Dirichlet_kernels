% --------------------- Pedro Saa UC 2021 ---------------------------------
addpath('data')

% Run 3 tuning methods with three different data sets and 3 components
clear
load settings1.mat
population11 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1); % orginal dataset
population21 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population31 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 5*xdata;                                                                        % increased data 5-fold, tol = 10%-total_data
population12 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population22 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population32 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 2*xdata;                                                                        % increased data 10-fold, tol = 10%-total_data
population13 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population23 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population33 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
save tunning_results_1

% Run 3 tuning methods with three different data sets and 4 components
clear
load settings2.mat
population11 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1); % orginal dataset
population21 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population31 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 5*xdata;                                                                        % increased data 5-fold, tol = 10%-total_data
population12 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population22 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population32 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 2*xdata;                                                                        % increased data 10-fold, tol = 10%-total_data
population13 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population23 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population33 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
save tunning_results_2

% Run 3 tuning methods with three different data sets and 5 components
clear
load settings3.mat
population11 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1); % orginal dataset
population21 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population31 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 5*xdata;                                                                        % increased data 5-fold, tol = 10%-total_data
population12 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population22 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population32 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 2*xdata;                                                                        % increased data 10-fold, tol = 10%-total_data
population13 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population23 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population33 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
save tunning_results_3

% Run 3 tuning methods with three different data sets and 6 components
clear

load settings4.mat
population11 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1); % orginal dataset
population21 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population31 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 5*xdata;                                                                        % increased data 5-fold, tol = 10%-total_data
population12 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population22 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population32 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
xdata = 2*xdata;                                                                        % increased data 10-fold, tol = 10%-total_data
population13 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),1);
population23 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),2);
population33 = multinomial_example(xdata,prior,N,alpha,tolStart,ceil(sum(xdata)/10),3);
save tunning_results_4