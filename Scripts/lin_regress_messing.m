x = linspace(0.001,1,300);
x = x';


%polynomial linear regression
X = [ones(length(x),1) x, x.^2, x.^3];
b = X\y ;
plot(x, X*b)
hold on 
plot(x,y)

%reciprocal linear regression
X = [ones(length(x),1) 1./x, 1./(x.^2)];
b = X\y ;

plot(1./x, 1/(X*b))
hold on 
plot(1./x,1/y)

figure()
plot(x,y)
hold on 
plot(x, 1./(X*b))



%log 
X = [ones(length(x),1) log(x), log(x.^2)];
b = X\y ;
figure()
plot(log(x), log(y))
plot(log(x), X*b)
