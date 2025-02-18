clear all
close all
clc

samplerate = 1.28;

%eps_values = [1,3,5,10,12,15].*10^(-4);
eps_values = 1e-3;
%eps = 1e-3;
%alpha_values = [0.0001, 0.001, 0.01, 0.1, 1, 10];
alpha = 0.01;
N_pts_per_wave=1;
N_windows=1; % Number of windows for comparison
N_MC=1000; % Number of Monte Carlo samples while solving the nonlinear function optimization for frequencies (4)
%num_modes = [4,6,8,12,16,20,24];
N_fmodes=16; % Number of elementary waves/ Fourier modes, i.e. parameter Nf in equation % (1)

Num_Windows = 1;

t = readmatrix('Time_DiabloCanyon_3_3.txt');
ele_val = readmatrix('DiabloCanyon_3_3.txt');
t_interval = [t(1),t(end)];
t_new = t(1):1/samplerate:t(end);
method = 'linear';
ele_new = interp1(t,ele_val,t_new, method);

figure;
f=gcf;
ax=gca;
plot(t, ele_val,'k-.','Linewidth',1.0)
hold on;
plot(t_new, ele_new,'r-','Linewidth',2.0)
yline(mean(ele_new),'b-.','Linewidth',1.5)
xlabel('Time (in seconds)','fontname','helvetica','fontsize',16);
ylabel('Surface elevation (\eta)','fontname','helvetica','fontsize',16);
xlim([min(min(t),min(t_new)) max(max(t),max(t_new))]);
ax.FontSize = 16;
name = sprintf('Original and interpolated waves.tiff');
legend('Original wave','Interpolated wave','Mean wave value','location','southeast','NumColumns',3,'fontSize',12);
grid on;
f.PaperUnits = 'inches';
f.PaperPosition = [0 0 15 4];
print(f,name,"-dtiff",'-r600'); 

% for i=1:length(eps_values)
% %for i=1:length(alpha_values)
% %for i=1:length(num_modes)
%     eps=eps_values(i);
%     %alpha = alpha_values(i);
%     %N_fmodes = num_modes(i);
%     
%     h_tmp=cell(1, 1);
%     t_tmp=cell(1, 1);
%     hz_slowflow=cell(1, 1);
%     slow_vars=cell(1, 1);
%     
%     hz = ele_new;
%     tz = 1:length(t_new);
%     
%     % extract extreme values, i.e. trough depths and crests heights
%     %[tmp_wave_height ,tmp_wave_idx, zero_idx] = my_wave_height_filter(hz ,N_pts_per_wave);
% %     [tmp_wave_height ,tmp_wave_idx] = my_wave_height_filter2(hz);
% %     
% %     ts=[tmp_wave_idx];%  tmp_t0];%[1:window_size];%
% %     hs=[tmp_wave_height];% zeros(1,length(zero_idx))];%hz(jj,:);%
%     
%     [ts,id]=sort(tz);
%     hs=hz(id);
%     
%     N_waves_tot = length(hs);
%     t_full=ts(1:N_waves_tot);
%     h_full=hs(1:N_waves_tot);
% 
%     t_full_plot = (t_full(1):t_full(end));
%     h_full_plot = hz(t_full(1):t_full(end));
% 
%     if N_fmodes<9
%         % If less than nine Fourier modes are selected a MATLAB routine can be
%         % used for fitting.
%         str_fmode=append('sin',num2str(N_fmodes));
%         my_fit=fit(t_full.',(h_full-mean(h_full)).',str_fmode);
%         err_const=sum((h_full-mean(h_full)-my_fit(t_full).').^2);
%         coeff2=coeffvalues(my_fit);
%         ws=sort(coeff2(2:3:end));
%         tic
%         % Monte Carlo sampling for optimization (4)
%         [ws_MC, ~]= my_Freq_MC_fit(t_full,h_full,err_const,N_fmodes,ws,N_MC,samplerate);
%         toc
%     else
%         ws=zeros(1,N_fmodes);
%         err_const=sum((h_full-mean(h_full)).^2);
%         tic
%         % Monte Carlo sampling for optimization (4)
%         [ws_MC, ~]= my_Freq_MC_fit2(t_full,h_full,err_const,N_fmodes,ws,N_MC,samplerate);
%         toc
%     end
%     
%     % Fitting the slowly varying amplitudes, i.e. solving the linear
%     % optimization (5), for the whole duration (cf. Fig. 3).
%     tic
%     h_tmp{1}= h_full;
%     t_tmp{1}= t_full;
%     [hz_slowflow{1} ,slow_vars{1}, ~]= my_SF_fit(t_full, h_full,ws_MC,eps,alpha);
%     toc
%     
%     % Utilize the slowly varying amplitudes from the whole time interval to
%     % fit the dynamical systems to be used for comparisons (cf. Fig. 3).
%     slow_vars_fit=slow_vars{1};
%     tt=t_tmp{1}(1):t_tmp{1}(end);
%     
%     % Interpolation yielding uniform time series for the slowly varying amplitudes.
%     slow_vars_interp_full = interp1(t_tmp{1},slow_vars_fit.',tt,'spline','extrap').';
%     writematrix(slow_vars_interp_full, 'SlowVaryingAmplitudes.csv');
% 
%     %name_amp = sprintf('Slow_amp_whole_%i.csv',i);
%     %name_time = sprintf('Whole_Time_%i.csv',i);
%     %writematrix(slow_vars_interp_full,name_amp,'Delimiter',',');
%     %writematrix(t_tmp{1}(1):t_tmp{1}(end),name_time,'Delimiter',',');
% %     hz_slowflow_full_truth=SF2fulltime(t_tmp{1},t_tmp{1},slow_vars{1},ws_MC,N_fmodes);
%     
% %     name_fit = sprintf('Fitted_wave_alpha=%d.txt',alpha);
% %     name_fit_time = sprintf('Fitted_wave_time_alpha=%d.txt',alpha);
% %     name_wave = sprintf('Original_wave_alpha=%d.txt',alpha);
% %     name_wave_time = sprintf('Original_wave_time_alpha=%d.txt',alpha);
%     
% %     name_fit = sprintf('Fitted_wave_Nf=%d.txt',N_fmodes);
% %     name_fit_time = sprintf('Fitted_wave_time_Nf=%d.txt',N_fmodes);
% %     name_wave = sprintf('Original_wave_Nf=%d.txt',N_fmodes);
% %     name_wave_time = sprintf('Original_wave_time_Nf=%d.txt',N_fmodes);
%     
% %     name_fit = sprintf('Fitted_wave_eps=%d.txt',eps);
% %     name_fit_time = sprintf('Fitted_wave_time_eps=%d.txt',eps);
% %     name_wave = sprintf('Original_wave_eps=%d.txt',eps);
% %     name_wave_time = sprintf('Original_wave_time_eps=%d.txt',eps);
% % 
% %     %dlmwrite(name_fit,hz_slowflow_full_truth);
% %     %dlmwrite(name_fit_time,t_full/(60 * samplerate));
% %     %dlmwrite(name_wave, h_full_plot);
% %     %dlmwrite(name_wave_time,t_full_plot/(60 * samplerate));
% % 
% %     figure;
% %     f=gcf;
% %     ax=gca;
% % 
% %     %ax.FontName = 'Helvetica';
% %     plot(t_full/(60 * samplerate), hz_slowflow_full_truth,'k','Linewidth',1.0)
% %     hold on;
% %     plot(t_full_plot/(60 * samplerate), h_full_plot,'r-','Linewidth',0.5)
% % 
% %     xlabel('Time (in seconds)','fontname','helvetica','fontsize',16);
% %     ylabel('Surface elevation (\eta)','fontname','helvetica','fontsize',16);
% % 
% %     xlim([t_full_plot(1)/(60 * samplerate) t_full_plot(end)/(60 * samplerate)]);
% %     ax.FontSize = 16;
% %     %name = sprintf('Wave_fit_model_predictions_alpha=%d.tiff',alpha);
% %     %name = sprintf('Wave_fit_model_predictions_Nf=%d.tiff',N_fmodes);
% %     name = sprintf('Wave_fit_model_predictions_eps=%d.tiff',eps);
% %     legend('Fit','real wave data','location','south','NumColumns',2,'fontSize',12);
% %     %name2 = sprintf('Wave fit alpha=%d',alpha);
% %     %name2 = sprintf('Wave fit Nf=%d',N_fmodes);
% %     name2 = sprintf('Wave fit eps=%d',eps);
% %     title(name2,'fontSize',12);
% %     grid on;
% %     f.PaperUnits = 'inches';
% %     f.PaperPosition = [0 0 15 4];
%     %print(f,name,"-dtiff",'-r600');
%     %print(f,name,"-dtiff",'-r600');
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% %pyrunfile('ML_Models.py')