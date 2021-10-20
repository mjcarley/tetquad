z="0.01" ;
z="1" ;
##z="0.1" ;

bdata = load(["basic-error-" z ".dat"]) ;
adata = load(["adaptive-error-" z ".dat"]) ;
ddata = load(["duffy-error-" z ".dat"]) ;

nb = bdata(:,1) ; eb = bdata(:,2) ;
nd = ddata(:,1) ; ed = ddata(:,2) ;

d = adata(:,3) ;
na = adata(:,1).*(4.^d-1)/3 ;
ea = adata(:,2) ;
[na,i] = sort(na) ; ea = ea(i) ;
