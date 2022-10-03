function Aq=interp2fun(x,y,A,xq,yq)

[xm,ym]=meshgrid(x,y);
[xqm,yqm]=meshgrid(xq,yq);
Aq=interp2(xm,ym,A,xqm,yqm,'spline');



end


