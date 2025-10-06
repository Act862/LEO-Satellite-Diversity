x = linspace(-3,3,200);
y = sin(x) + 0.1*randn(size(x));

x_train = x(1:2:end);
y_train = y(2:2:end);
x_test = x(2:2:end);
y_test = y(2:2:end);

net = feedforwardnet(10);
net = train(net,x_train,y_train);

y_pred = net(x_test);

figure;
plot(x_train,y_train,'bo',x_test,y_test,'r-');
legend('Training Data', 'Neural Network Output');
title('Function Approximation');