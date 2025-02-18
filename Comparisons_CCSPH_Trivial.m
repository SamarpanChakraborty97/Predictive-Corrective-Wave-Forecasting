clc
clear all

elevation = readmatrix("DiabloCanyon_3_3.txt");
time = readmatrix("Time_DiabloCanyon_3_3.txt");

num_gauges = length(elevation(1,:));
num_entries = length(elevation(:,1));

Gauges = 3;

error_c = zeros(length(Gauges),1);
error_n = zeros(length(Gauges),1);

for j=1:num_gauges
       
    t = time(:,j);
    
    mean_c = mean(elevation(:,j));
    ele_c = elevation(:,j) - mean_c;

    f = fit(t,ele_c,'smoothingspline','SmoothingParam',0.99);
    h = plot(f);
    xi_c = get(h,'XData');
    yi_c = get(h,'YData');
    hold off;
 
    figure;
    g = gcf;
    ax = gca;
    ax.FontSize = 15;
    plot(xi_c,yi_c,'k-.','Linewidth',2)
    grid on;
    hold off;
    
    name = sprintf("Displacement at gauge location");
    title(name, 'FontName','Times','FontSize',13, 'Fontweight','normal');
    xlabel('Time (seconds)','FontSize', 12)
    ylabel('\eta (metres)','FontSize',12)
    
    g.PaperUnits = 'inches';
    g.PaperPosition = [0 0 16.5 3];
    name2 = sprintf('Model predictions at gauge location_3.png');
    saveas(g,name2)
    
    close all;
end
