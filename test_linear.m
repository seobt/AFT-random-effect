function fun = test_linear (b,x)
p=length(b);
fun=b(1)+x(:,1:p-1)*b(2:p);
