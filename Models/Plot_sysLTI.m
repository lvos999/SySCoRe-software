function Plot_sysLTI(sysLTI)
%PLOT_SYSLTI Plot the important regions of the LTI system

if size(sysLTI.C,1)~=2
    error('This function does not work for systems with a non-2D space')
end
%1. plot the domain of X
if sysLTI.dim==2
    figure;
    plot_x = plot(sysLTI.X);
    set(plot_x, 'FaceColor','None');
    title('X domain')
end

%2. plot the regions of the LTI system together with the AP strings
figure;
plot_x = plot(sysLTI.C*sysLTI.X);
set(plot_x, 'FaceColor','None');
hold on 
for i =1:length(sysLTI.regions)
    plot_x = plot(sysLTI.regions(i));
    set(plot_x, 'FaceColor','None');
    hold on
    xc = sysLTI.regions(i).chebyCenter;
    text(xc.x(1),xc.x(2), sysLTI.AP{i})
end

title('Output space Y with labelled regions' )



end

