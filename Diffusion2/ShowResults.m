clear
close all

% load 4DVarTest
% load ETKFTest

ArelRMSE_bg = averageRelativeRootMeanSquareError(u_bg_tot, u_tot)
ArelRMSE_meas = averageRelativeRootMeanSquareError(u_meas, u_tot(2:nx-1,2:end))
ArelRMSE_a = averageRelativeRootMeanSquareError(u_a_tot, u_tot)

f = figure;

for i=1    
    subplot(3,1,1)
    plot(x, u_tot(:,i));
    suptitle({['1-D Diffusion with \nu =',num2str(vis),' and \beta = ',num2str(beta)];['time(\itt) = ',num2str(dt*(i-1))]})
    axis([0 2 0 3])
    
    subplot(3,1,2)
    plot(x, u_a_tot(:,i));
    axis([0 2 0 3])
    ylabel('Transport property profile (u) \rightarrow')
    
	subplot(3,1,3)    
	plot(x, u_bg_tot(:,i));    
    axis([0 2 0 3])
    xlabel('Spatial co-ordinate (x) \rightarrow')
    drawnow; 
    refreshdata(f)
end
pause(0.3);
for i=2:nt+1
    subplot(3,1,1)
    suptitle({['1-D Diffusion with \nu =',num2str(vis),' and \beta = ',num2str(beta)];['time(\itt) = ',num2str(dt*(i-1))]})
    plot(x, u_tot(:,i), x(2:end-1), u_meas(:, i-1), '*');
    axis([0 2 0 3])    	
    
    subplot(3,1,2)
    plot(x, u_a_tot(:,i));
    axis([0 2 0 3])
    ylabel('Transport property profile (u) \rightarrow')
    
	subplot(3,1,3)    
	plot(x, u_bg_tot(:,i));    
    axis([0 2 0 3])
    xlabel('Spatial co-ordinate (x) \rightarrow')
    
    drawnow; 
    refreshdata(f)
    pause(dt);
end