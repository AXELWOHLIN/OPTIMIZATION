disp("hi")


%%%%% Axel Wohlin Optimization 2 MatLab code

expected_returns = 0.01*[13; 4.4; 12.1; 7.1; 11.7];

Sigma = [
    4.01, -1.19, 0.60, 0.74, -0.21;
    -1.19, 1.12, 0.21, -0.54, 0.55;
    0.6, 0.21, 3.31, 0.77, 0.29;
    0.74, -0.54, 0.77, 3.74, -1.04;
    -0.21, 0.55, 0.29, -1.04, 2.6] * 0.01;



%%%%%%%%%%%%%%%% Task 2
A11 = Sigma;
% Number of assets
n = 5;
A12 = -ones(n,1);
A21 = -A12';
A22 = 0;

b1 = zeros(n,1);
b2 = 1;

A = [A11,A12;
    A21,A22];

b = [b1;b2];

sol = A\b;
w = sol(1:size(Sigma, 1));
disp("Task 2!")
variance = w' * Sigma * w;
disp("Optimal weights:");
disp(w);
disp("Variance:");
disp(variance);
disp("Expected ROI")
disp(w'*expected_returns)


%%%%%%%% Task 3 no short selling
disp("Task 3 without short sell")
return_value = 0.12;

A11 = Sigma;
% Number of assets
n = 5;
A12 = -[ones(n,1),expected_returns];
A21 = -A12';
A22 = zeros(2,2);

b1 = [zeros(n,1),ones(n,1)];
b2 = [1,return_value];

A = [A11,A12;
    A21,A22];

b = [b1;b2];

sol = quadprog(A11, [], [], [], A21, b2, zeros(n,1) );
w = sol(1:size(Sigma, 1));
variance = w' * Sigma * w;
disp("Optimal weights:");
disp(w);
disp("Variance:");
disp(variance);
disp("Expected ROI")
disp(w'*expected_returns)


%%%%%%%% Task 3 with short selling
disp("Task 3 with short sell")
return_value = 0.12;

A11 = Sigma;
% Number of assets
n = 5;
A12 = -[ones(n,1),expected_returns];
A21 = -A12';
A22 = zeros(2,2);

b1 = [];
b2 = [1,return_value];

A = [A11,A12;
    A21,A22];

b = [b1;b2];

sol = quadprog(A11, [], [], [], A21, b2, []);
w = sol(1:size(Sigma, 1));
variance = w' * Sigma * w;
disp("Optimal weights:");
disp(w);
disp("Variance:");
disp(variance);
disp("Expected ROI")
disp(w'*expected_returns)


%%%%% Task 4 with short sell
disp("Task 4 with short sell!")


returns = [];
variances = [];
A11 = Sigma;
A12 = -[ones(size(Sigma, 1),1)];
A21=-A_12';
A22=0;
b1 = zeros(size(Sigma, 1), 1);
b2 = 1;
A = [A11, A12; A21, A22];
b = [b1; b2];
for alpha=0.05:0.05:1
    frontier = -1*(1-alpha)*expected_returns;
    sol = quadprog(A_11, frontier, [], [], A21, b2, [] );
    w = sol(1:size(Sigma, 1));
    return_value= w'*expected_returns;
    variance = w' * Sigma * w;
    returns=[returns,return_value];
    variances=[variances,variance];
end
hold on
figure(1)
plot(variances, returns*100)
xlabel('Variance')
ylabel('Expected ROI [%]')
title("Plot over alpha with short sell")


%%%%Task 4 without short sell
disp("Task 4 without short sell")

returns = [];
variances = [];
A11 = Sigma;
A12 = -[ones(size(Sigma, 1),1)];
A21=-A_12';
A22=0;
b1 = zeros(size(Sigma, 1), 1);
b2 = 1;
A = [A11, A12; A21, A22];
b = [b1; b2];
for alpha=0.05:0.05:1
    frontier = -1*(1-alpha)*expected_returns;
    sol = quadprog(A_11, frontier, [], [], A21, b2, zeros(n,1) );
    w = sol(1:size(Sigma, 1));
    return_value= w'*expected_returns;
    variance = w' * Sigma * w;
    returns=[returns,return_value];
    variances=[variances,variance];
end
hold on
figure(2)
plot(variances, returns*100)
xlabel('Variance')
ylabel('Expected ROI [%]')
title("Plot over alpha without short sell")