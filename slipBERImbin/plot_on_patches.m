% Script to plot which patches on a fault plane are 'on' and which are
% 'off', plus the circular harmonics shape that's defining which are on/off
%
% rmja 23-apr-2018

circharm_center_trial = m_curr(m_identifyer_master==9);

figure('position', [400, 800, 1200, 1000]);
scatter(patchx(onoffidentifyer==1), patchz(onoffidentifyer==1), 1800, 'gs', 'filled')   % plot x and z of on patches
hold on;
scatter(patchx(onoffidentifyer==0), patchz(onoffidentifyer==0),1800, 'rs', 'filled')   % plot x and z of off patches
plot(circx+circharm_center_trial(1),circz+circharm_center_trial(2));
set(gca,'Ydir','reverse')
axis equal tight
axis([min([circz+circharm_center_initial(1) 0]) max([circz+circharm_center_initial(1)  fault_length_for_smoothing]) min([circx+circharm_center_initial(2) 0]) max([circx+circharm_center_initial(2)  fault_width_for_smoothing])])
ylabel('Down dip')
xlabel('Along strike')
title('Slipping area on fault')