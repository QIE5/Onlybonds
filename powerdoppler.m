
%SVD parameter
nCut = 4;
nNoiseCut = 0;

%SVD clutter filtering
data_filt = svdClutterFilter(data, nCut, nNoiseCut);

%Power Doppler
PD = sum(abs(data_filt).^2, 3);

%convert to dB
PD_dB = 10*log10(PD + eps);
PD_dB = PD_dB - max(PD_dB(:));


PD_dB = rot90(PD_dB, 2);


pixelSize = 3.3879e-5;


[rows, cols] = size(PD_dB);

dx = pixelSize;
dz = pixelSize;


xAxis_mm = ((1:cols) - 0.5) * dx * 1e3;
zAxis_mm = ((1:rows) - 0.5) * dz * 1e3;


figure;
imagesc(xAxis_mm, zAxis_mm, PD_dB, [-40 0]);
set(gca, 'YDir', 'normal');
colormap hot;
colorbar;
axis image off;

figure;
imagesc(xAxis_mm, zAxis_mm, PD_dB, [-40 0]);

colormap hot;
colorbar;
axis image off;


set(gca, 'YDir', 'normal');

hold on;


barLength_mm = 1;
marginFrac   = 0.06;

ax = gca;
xl = xlim(ax);
yl = ylim(ax);

xRange = xl(2) - xl(1);
yRange = yl(2) - yl(1);


x0 = xl(1) + marginFrac * xRange;


if strcmp(ax.YDir, 'normal')
    y0 = yl(1) + marginFrac * yRange;
else
    y0 = yl(2) - marginFrac * yRange;
end


plot([x0, x0 + barLength_mm], [y0, y0], ...
    'w', 'LineWidth', 4, 'Clipping', 'on');

hold off;

exportgraphics(gcf, 'power_doppler.png', 'Resolution', 300);