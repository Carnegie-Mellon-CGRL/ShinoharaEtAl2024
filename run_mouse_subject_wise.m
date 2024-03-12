clear all; clc; clf;

global connectivity diameter

connectivity=[...
    0.19 4.06 2.48 1.36 0 0 0 0 0 0 0;
    0 0.12 1.55 0.93 0 0 0 0 0 0 0;
    0 0 0.17 2.26 0.31 0.05 0.11 0.13 0 0 0;
    0 0 0 0.05 2.00 0.70 0.62 0.43 0.39 0 0;
    0 0 0 0 0.18 1.92 1.56 0.85 0.83 0.50 0;
    0 0 0 0 0 0.18 2.05 1.15 1.00 0.67 0;
    0 0 0 0 0 0 0.05 1.55 1.06 0.67 0;
    0 0 0 0 0 0 0 0.08 1.67 0.83 6;
    0 0 0 0 0 0 0 0 0.06 1.67 4;
    0 0 0 0 0 0 0 0 0 0.17 2;
    0 0 0 0 0 0 0 0 0 0 0;
    ];

connectivity = connectivity(1:10, 1:10);

%Initialize vessel order metrics
n_orders = size(connectivity,1); %Number of vessel orders, vessels decrease in size with decreasing order
orders = 1:n_orders;

filein = 'input_data.txt';
delimiter = '\t';
input_data = importdata(filein, delimiter);

%Parametric arrays
%[Sham Shunt ShuntTX]
exp_groups = ["Sh"; "Smt"; "ACS"; "Tx"];
n_groups = size(exp_groups,1);
per_group = zeros(n_groups,1);
n_samples = size(input_data.textdata,1);
sample_group_type = zeros(n_samples,1);
for i = 1:n_samples
    for j = 1:n_groups 
        has_group_name = strfind(input_data.textdata{i}, exp_groups{j});
        if size(has_group_name,1) > 0
            sample_group_type(i) = j;
            per_group(j) = per_group(j) + 1;
        end
    end
end

r_iv_LPA_sample = input_data.data(:,1) * 10^-6; %m
diam_ratio_const_sample = input_data.data(:,2); %ratio
CO_sample = input_data.data(:,3) / 60; %mL / s
SPAP_sample = input_data.data(:,4); %mmHg
DPAP_sample = input_data.data(:,5); %mmHg
LVEDP_sample = input_data.data(:,6); %mmHg
FlowSplit_sample = input_data.data(:,7); %ratio

group_count = zeros(n_groups,1);
hemo_groups = cell(n_groups,1);

for i = 1:n_samples
    display("Sample: " + int2str(i))
    input_pars.n_orders = n_orders;
    input_pars.r_iv_LPA = r_iv_LPA_sample(i);
    input_pars.h_iv_LPA = 40 * 10^-6;
    input_pars.l_iv_LPA = 0.0053;
    input_pars.V_LPA = 4.39 * 10^-10;
    input_pars.len_ratio = 1.60;
    
    %Average measurements for sham animals
    SPAP = SPAP_sample(i); %mmHg
    DPAP = DPAP_sample(i); %mmHg
    LVEDP = LVEDP_sample(i); %mmHg
    CO = CO_sample(i); %mL / s

    %Derived
    Q_LPA = CO * FlowSplit_sample(i);
    P_term = LVEDP;
    mPAP = (2 / 3 * DPAP + 1 / 3 * SPAP) * 1333.22; %dyne/cm^2
    dP = mPAP - LVEDP * 1333.22; %dyne/cm^2
    PVR = dP / CO; %cgs resistance
    PVR_LPA = dP / Q_LPA;
    input_pars.PVR_LPA = PVR_LPA;
    
    lb = 1.0;
    ub = 2.0;
    par0 = diam_ratio_const_sample(i);
    options = optimset('MaxFunEvals',16000,'MaxIter',16000,'TolFun',1e-12,'TolX',1e-12,...
                           'Display','iter');
    current_objective = @(par)ord_diam_obj(par, input_pars);
    par_est = lsqnonlin(current_objective, par0, lb, ub, options);
    diam_ratio_const_sample(i) = par_est;
    
    %Calculate hemodynamics for given orders
    order = n_orders;
    order_remodeling = n_orders;
    alpha = 1;
    beta = 1;
    PVR_LPA_tree = find_opt_morphometric_tree(order, order_remodeling, alpha, beta);

    hemo = calculate_hemo_morphometric_tree_termBC(Q_LPA, P_term);
    results.(input_data.textdata{i}).diameter = diameter;
    results.(input_data.textdata{i}).hemo = hemo;
    results.(input_data.textdata{i}).PVR_LPA = PVR_LPA;
    results.(input_data.textdata{i}).PVR_LPA_tree = PVR_LPA_tree;
    
    group_count(sample_group_type(i)) = group_count(sample_group_type(i)) + 1;
    
    %Store hemodynamic data from each group for stats
    hemo_groups{sample_group_type(i)}.P_order_all(: , group_count(sample_group_type(i))) = hemo.P_order(:,1);
    hemo_groups{sample_group_type(i)}.Q_order_all(: , group_count(sample_group_type(i))) = hemo.Q_order(:,1);
    P_order_outlet = [P_term * 1333.22; hemo.P_order(1:end - 1,1)];
    hemo_groups{sample_group_type(i)}.Res_order_all(:, group_count(sample_group_type(i))) = (hemo.P_order(:,1) - P_order_outlet) ./ hemo.Q_order(:,1);
    hemo_groups{sample_group_type(i)}.Sigma_order_all(: , group_count(sample_group_type(i))) = hemo.Sigma_order(:,1);
    hemo_groups{sample_group_type(i)}.WSS_order_all(: , group_count(sample_group_type(i))) = hemo.WSS_order(:,1);
    hemo_groups{sample_group_type(i)}.diameter_all(: , group_count(sample_group_type(i))) = diameter;
    
%     fig_num = i;
%     plot_hemo(hemo, order, fig_num)
end

%Calculate means for each group
for j = 1:n_groups
    %averages
    hemo_groups{j}.P_order(:,1) = mean(hemo_groups{j}.P_order_all,2);
    hemo_groups{j}.Q_order(:,1) = mean(hemo_groups{j}.Q_order_all,2);
    hemo_groups{j}.Res_order(:,1) = mean(hemo_groups{j}.Res_order_all,2);
    hemo_groups{j}.Sigma_order(:,1) = mean(hemo_groups{j}.Sigma_order_all,2);
    hemo_groups{j}.WSS_order(:,1) = mean(hemo_groups{j}.WSS_order_all,2) ;
    hemo_groups{j}.diameter(:,1) = mean(hemo_groups{j}.diameter_all,2);
    %stds
    hemo_groups{j}.P_order(:,2) = std(hemo_groups{j}.P_order_all,0,2);
    hemo_groups{j}.Q_order(:,2) = std(hemo_groups{j}.Q_order_all,0,2);
    hemo_groups{j}.Res_order(:,2) = std(hemo_groups{j}.Res_order_all,0,2);
    hemo_groups{j}.Sigma_order(:,2) = std(hemo_groups{j}.Sigma_order_all,0,2);
    hemo_groups{j}.WSS_order(:,2) = std(hemo_groups{j}.WSS_order_all,0,2);
    hemo_groups{j}.diameter(:,2) = std(hemo_groups{j}.diameter_all,0,2);

end

hemo_mult = [hemo_groups{1} hemo_groups{2} hemo_groups{3} hemo_groups{4}];

fig_num = i + 1;
plot_hemo_mult(hemo_mult, order, fig_num)

colors = [0 1 0
          0.8 1 0.8
          1 0 0
          1 0.8 0.8];

fig_num = i + 2;
figure(fig_num)
for i = 1:3
     
     sp2 = semilogy(1:n_orders, hemo_groups{i}.diameter(:,1) * 10^4, 'LineWidth',1.5, 'Color', colors(i,:));
     xlabel('Vessel Order (-)'); ylabel('Diameter (um)')
     legend('Sham', 'Sham+Tx', 'Shunt', 'Shunt+Tx')
     hold on
     set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
    
end

save('morphometric_tree_results.mat', 'results', 'hemo_groups')





