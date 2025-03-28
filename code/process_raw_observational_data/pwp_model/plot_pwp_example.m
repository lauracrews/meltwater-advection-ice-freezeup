%
%
%
%




load example_Output

% Check heat conservation:
figure(1)
clf
plot(time,rho.*cp.*(sum(temp).*dz-sum(temp(:,1)).*dz))
hold on
plot(time,cumsum(rhf(tstart:tstop)+thf(tstart:tstop)).*dt,'r')
title('Heat content')
legend('Water column heat content', 'Surface heat input')
ylabel('J/m^2')
gregaxd(time,1)


%%

figure(2)
clf
subplot(221)
contourf(time,z,temp,[27.5:.25:30.],'linestyle','none')
set(gca,'ydir','reverse')
ax=axis;
axis([ax(1:2) 0 20])
ax=axis;
caxis([27.5 30])
title('Model temperature')
gregaxd(time,1)



subplot(222)
contourf(data.time,data.z,ztemp1n',[27.5:.25:30.],'linestyle','none')
set(gca,'ydir','reverse')
axis(ax)
caxis([27.5 30])
title('Observed temperature')
gregaxd(time,1)
colorbar_out;

subplot(223)
contourf(time,z,temp,[24.5:.25:30.],'linestyle','none')
set(gca,'ydir','reverse')
ax=axis;
axis([ax(1:2) 0 45])
ax=axis;
caxis([24.5 30])
gregaxd(time,1)
title('Model temperature')

subplot(224)
contourf(ztimetemp1n,zdepthtemp1n,ztemp1n',[24.5:.25:30.],'linestyle','none')
set(gca,'ydir','reverse')
axis(ax)
caxis([24.5 30])
gregaxd(time,1)
colorbar_out;
title('Observed temperature')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
znorth_low=run_avg2(znorth1n,13./24,1);
zeast_low=run_avg2(zeast1n,13./24,1);
znorth_bp=run_avg2(znorth1n,13./24,1)-run_avg2(znorth1n,13,1);
zeast_bp=run_avg2(zeast1n,13./24,1)-run_avg2(zeast1n,13,1);
znorth_low_rel=znorth_low-znorth_low(:,7)*ones(1,length(zdepthvel1n));
zeast_low_rel=zeast_low-zeast_low(:,7)*ones(1,length(zdepthvel1n));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
clf
subplot(221)
contourf(time,z,U,[-.5:.02:.5],'linestyle','none')
set(gca,'ydir','reverse')
ax=axis;
axis([ax(1:2) 0 45])
ax=axis;
caxis([-.5:.5])
gregaxd(time,1)
title('Model eastward velocity')

subplot(222)
contourf(time,z,U-ones(size(z))*U(21,:),[-.2:.02:.2],'linestyle','none');
set(gca,'ydir','reverse')
ax=axis;
axis([ax(1:2) 0 10])
ax=axis;
caxis([-.2 .2])
gregaxd(time,1)
title('Model eastward velocity, relative to 10.25 m')
colorbar_out;

subplot(223)
contourf(ztimevel1n,zdepthvel1n,zeast_low'./100,[-.5:.02:.5],'linestyle','none');
set(gca,'ydir','reverse')
axis(ax)
caxis([-.5:.5])
gregaxd(time,1)
title('Observed eastward velocity (smoothed)')

subplot(224)
contourf(ztimevel1n,zdepthvel1n,zeast_low_rel'./100,[-.5:.02:.5],'linestyle','none');
set(gca,'ydir','reverse')
axis(ax)
caxis([-.5:.5])
gregaxd(time,1)
title('Observed eastward velocity relative to 50 m')
colorbar_out;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
clf
subplot(221)
contourf(time,z,V,[-.5:.02:.5])
shading flat
set(gca,'ydir','reverse')
ax=axis;
axis([ax(1:2) 0 45])
ax=axis;
caxis([-.5:.5])
gregaxd(time,1)
title('Model northward velocity')

subplot(222)
contourf(time,z,V-ones(size(z))*V(21,:),[-.2:.02:.2])
shading flat
set(gca,'ydir','reverse')
ax=axis;
axis([ax(1:2) 0 10])
ax=axis;
caxis([-.2 .2])
gregaxd(time,1)
title('Model northward velocity, relative to 10.25 m')

subplot(223)
contourf(ztimevel1n,zdepthvel1n,znorth_low'./100,[-.5:.02:.5])
shading flat
set(gca,'ydir','reverse')
axis(ax)
caxis([-.5:.5])
gregaxd(time,1)
title('Observed northward velocity (smoothed)')

subplot(224)
contourf(ztimevel1n,zdepthvel1n,znorth_low_rel'./100,[-.5:.02:.5])
shading flat
set(gca,'ydir','reverse')
axis(ax)
caxis([-.5:.5])
gregaxd(time,1)
title('Observed northward velocity relative to 50 m depth (smoothed)')




