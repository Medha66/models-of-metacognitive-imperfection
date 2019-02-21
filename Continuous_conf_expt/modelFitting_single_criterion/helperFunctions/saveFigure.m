function [  ] = saveFigure(FigureName)

set(gcf,'PaperUnits','inches','PaperSize',[8 12],'PaperPosition',[1 1 6.65 5])
print(gcf,FigureName,'-djpeg')
delete(gcf)
end

