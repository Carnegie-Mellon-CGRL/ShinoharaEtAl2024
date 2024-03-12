function fig = plot_hemo_mult(hemo_all, n_orders, fig_num)

 [temp, num_hemo] = size(hemo_all);
 n_orders_all = ones(n_orders, num_hemo);
 n_orders_all = n_orders_all .* [1:n_orders]';
 
 P_all = zeros(n_orders, num_hemo);
 P_err_all = zeros(n_orders, num_hemo);
 Q_all = zeros(n_orders, num_hemo);
 Q_err_all = zeros(n_orders, num_hemo);
 Res_all = zeros(n_orders, num_hemo);
 Res_err_all = zeros(n_orders, num_hemo);
 Sigma_all = zeros(n_orders, num_hemo);
 WSS_all = zeros(n_orders, num_hemo);
 WSS_err_all = zeros(n_orders, num_hemo);
 for i = 1:num_hemo
     P_all(:,i) = hemo_all(i).P_order(:,1) / (10 * 133.33);
     P_err_all(:,i) = hemo_all(i).P_order(:,2) / (10 * 133.33);
     Q_all(:,i) = hemo_all(i).Q_order(:,1) * 60;
     Q_err_all(:,i) = hemo_all(i).Q_order(:,2) * 60;
     Res_all(:,i) = hemo_all(i).Res_order(:,1) / (10 * 133.33) / 60;
     Res_err_all(:,i) = hemo_all(i).Res_order(:,2) / (10 * 133.33) / 60;
     Sigma_all(:,i) = hemo_all(i).Sigma_order(:,1) / 10000;
     WSS_all(:,i) = hemo_all(i).WSS_order(:,1);
     WSS_err_all(:,i) = hemo_all(i).WSS_order(:,2);
 end
 
 fig = figure(fig_num);
 subplot(2,2,1)
 hold on
%  er1 = errorbar(n_orders_all, P_all, P_err_all * 0, P_err_all);
%  er1(1).Color = [0 0 0];
%  er1(2).Color = [0 0 0];
%  er1(3).Color = [0 0 0];
%  er1(1).LineStyle = 'none';
%  er1(2).LineStyle = 'none';
%  er1(3).LineStyle = 'none';
 
 sp1 = bar(n_orders_all, P_all, 'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'LineWidth',0.5);
 sp1(1).FaceColor = [0 1 0];
 sp1(2).FaceColor = [0.8 1.0 0.8];
 sp1(3).FaceColor = [1 0 0];
 sp1(4).FaceColor = [1 0.8 0.8];
 xlabel('Vessel Order (-)'); ylabel('Pressure (mmHg)')

 ylim([0 40])
 set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 
 subplot(2,2,2)
 hold on
%  er2 = errorbar(n_orders_all, Q_all, Q_err_all * 0, Q_err_all);
%  er2(1).Color = [0 0 0];
%  er2(2).Color = [0 0 0];
%  er2(3).Color = [0 0 0];
%  er2(1).LineStyle = 'none';
%  er2(2).LineStyle = 'none';
%  er2(3).LineStyle = 'none';
 
 sp2 = bar(n_orders_all, Q_all, 'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'LineWidth',0.5);
 sp2(1).FaceColor = [0 1 0];
 sp2(2).FaceColor = [0.8 1.0 0.8];
 sp2(3).FaceColor = [1 0 0];
 sp2(4).FaceColor = [1 0.8 0.8];
 xlabel('Vessel Order (-)'); ylabel('Flow (mL/min)')

 set(gca, 'YScale', 'log', 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 
%  subplot(2,2,3)
%  sp3 = bar(n_orders_all, Sigma_all,'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'LineWidth',0.5);
%  xlabel('Vessel Order (-)'); ylabel('Circ. Stress (kPa)')
%  hold on
% %  er = errorbar(1:n_orders, hemo.WSS_order(:,1), 0 * hemo.WSS_order(:,2), hemo.WSS_order(:,2));
% %  er.Color = [0 0 0];
% %  er.LineStyle = 'none';
% %  ylim([0 20])
%  set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 subplot(2,2,3)
 sp3 = bar(n_orders_all, Res_all, 'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'LineWidth',0.5);
 sp3(1).FaceColor = [0 1 0];
 sp3(2).FaceColor = [0.8 1.0 0.8];
 sp3(3).FaceColor = [1 0 0];
 sp3(4).FaceColor = [1 0.8 0.8];
 xlabel('Vessel Order (-)'); ylabel('Resistance (mmHg*min/ml)')
%  ylim([0 40])
 set(gca, 'YScale', 'log', 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 
 subplot(2,2,4)
 hold on
%  er4 = errorbar(n_orders_all, WSS_all, WSS_err_all * 0, WSS_err_all);
%  er4(1).Color = [0 0 0];
%  er4(2).Color = [0 0 0];
%  er4(3).Color = [0 0 0];
%  er4(1).LineStyle = 'none';
%  er4(2).LineStyle = 'none';
%  er4(3).LineStyle = 'none';
 
 sp4 = bar(n_orders_all, WSS_all, 'FaceColor',[0 0 0],'EdgeColor',[1 1 1],'LineWidth',0.5);
 sp4(1).FaceColor = [0 1 0];
 sp4(2).FaceColor = [0.8 1.0 0.8];
 sp4(3).FaceColor = [1 0 0];
 sp4(4).FaceColor = [1 0.8 0.8];
 xlabel('Vessel Order (-)'); ylabel('Wall Shear Stress (dyne/cm^2)')
 set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
 
 set(gcf,'Units','inches','Position',[.25 .25 16 4],'PaperUnits','inches','PaperPosition',[.25 .25 16 4])

end