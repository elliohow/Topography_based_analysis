function pRFDisplay(outputStruct, x, y, z)
%% display pRF for single voxel

structX = outputStruct.x(:,1,1);
structY = outputStruct.y(1,:,1);

figure

subplot(2,2,1)
imagesc(structX, structY, permute(outputStruct.data{x, y, z}.pinv_pRF, [2,1]))
formatPlot()
title('pinv')
axis equal

subplot(2,2,2)
imagesc(structX, structY, permute(outputStruct.data{x, y, z}.Lasso_pRF, [2,1]))
formatPlot()
title('lasso')
axis equal

subplot(2,2,3)
imagesc(structX, structY, permute(outputStruct.data{x, y, z}.Ridge_pRF, [2,1]))
formatPlot()
title('ridge regression')
axis equal

subplot(2,2,4)
imagesc(structX, structY, permute(outputStruct.data{x, y, z}.svm_pRF, [2,1]))
formatPlot()
title('svm regression')
axis equal

end

%% helper function
function formatPlot()
% formatPlot - format plot and add v/h lines

hold('on')
line([0,0], get(gca,'ylim'),'color','w', 'linewidth',2);
line(get(gca,'xlim'),[0,0],'color','w', 'linewidth',2);
axis('image')
colormap(parula())
colorbar()
xlabel('Visual space (x)');
ylabel('Visual space (y)');

end
