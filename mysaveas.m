function mysaveas(fig,fn)
if isMatlab,
  saveas(fig,sprintf('%s.pdf',fn),'pdf');
  saveas(fig,sprintf('%s.jpg',fn),'jpg');
  saveas(fig,sprintf('%s.jpg',fn),'eps')
end;
