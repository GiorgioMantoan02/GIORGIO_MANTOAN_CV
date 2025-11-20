clear all
clc
close all
clear all

%% preparazione dati

filename = 'Gruppo14_Lab5.xlsx';
sheetname= 'Historical_Simulation';
data = readmatrix(filename, 'Sheet', sheetname);

% file: Gruppo14_Lab5.xlsx
% foglio: Historical_Simulation
data = data(5:end,2:4);

loss_apple = diff(data(:,1));
loss_amazon = diff(data(:,2));
loss_facebook = diff(data(:,3));
loss = [loss_apple, loss_amazon, loss_facebook];

mu = mean([loss_apple, loss_amazon, loss_facebook])';
sigma = std([loss_apple, loss_amazon, loss_facebook]);

V_1 = cov(loss_apple, loss_amazon); 
cov_1 = V_1(1,2); % covarianza apple - amazon
V_2 = cov(loss_apple, loss_facebook); 
cov_2 = V_2(1,2); % covarianza apple - facebook
V = cov(loss);


%% ptf1-2: APPLE - AMAZON, VaR = 28, p1 = 0.05 p2 = 0.01

VaR1 = 28;
p1 = 0.05;
p2 = 0.01;

w1 = soluzioni_VaR(mu(1), mu(2), sigma(1), sigma(2), cov_1, VaR1, p1);
w1_1 = [w1(1); 1-w1(1)];
w1_2 = [w1(2); 1-w1(2)];

w2 = soluzioni_VaR(mu(1), mu(2), sigma(1), sigma(2), cov_1, VaR1, p2);
w2_1 = [w2(1); 1-w2(1)];
w2_2 = [w2(2); 1-w2(2)];


%% ptf3-4: APPLE - FACEBOOK, VaR = 28, p3 = 0.05 p4 = 0.01

w3 = soluzioni_VaR(mu(1), mu(3), sigma(1), sigma(3), cov_2, VaR1, p1);
w3_1 = [w3(1);1-w3(1)];
w3_2 = [w3(2);1-w3(2)];

w4 = soluzioni_VaR(mu(1), mu(3), sigma(1), sigma(3), cov_2, VaR1, p2);
w4_1 = [w4(1);1-w4(1)];
w4_2 = [w4(2);1-w4(2)];


%% punto 5 : VaR 1% (t+1) su w_MVP

p = 0.01;

returns_apple = log(data(1:end-1,1)./data(2:end,1));
returns_amazon = log(data(1:end-1,2)./data(2:end,2));
returns_facebook = log(data(1:end-1,3)./data(2:end,3));
returns = [returns_apple, returns_amazon, returns_facebook];

V_ret = cov(returns);
Ones = ones(3,1);
C = Ones' * inv(V_ret) * Ones;

w_mvp = inv(V_ret) * Ones / C;

mu_mvp = w_mvp' * mu;
sigma_mvp = sqrt(w_mvp' * V * w_mvp);

VaR_mvp = mu_mvp + sigma_mvp * norminv(1-p);

ES_mvp = calcolo_ES(w_mvp, mu, V, p);

fprintf("\nw_MVP = [%.4f, %.4f, %.4f]\n", w_mvp)
fprintf("VaR_0.01(t+1) = %.4f\nES_0.01(t+1) = %.4f\n", VaR_mvp, ES_mvp)

%% expected shortfall
% ES di ogni w_i, solo per il rispettivo valore di p (usato per il calcolo
% del ptf tramite il VaR-target)

% ptf 1
ES_w1_1 = calcolo_ES(w1_1, mu([1 2]),V_1,p1);
ES_w1_2 = calcolo_ES(w1_2, mu([1 2]),V_1,p1);

% ptf 2
ES_w2_1 = calcolo_ES(w2_1, mu([1 2]),V_1,p2);
ES_w2_2 = calcolo_ES(w2_2, mu([1 2]),V_1,p2);

% ptf 3
ES_w3_1 = calcolo_ES(w3_1, mu([1 3]),V_2,p1);
ES_w3_2 = calcolo_ES(w3_2, mu([1 3]),V_2,p1);

% ptf 4
ES_w4_1 = calcolo_ES(w4_1, mu([1 3]),V_2,p2);
ES_w4_2 = calcolo_ES(w4_2, mu([1 3]),V_2,p2);

ES = [ES_w1_1, ES_w1_2; ES_w2_1, ES_w2_2; 
      ES_w3_1, ES_w3_2; ES_w4_1, ES_w4_2];

fprintf("\nPtf1 VaR = 28 p = %.2f:\n", p1)
fprintf("\tw1 = [%.4f, %.4f], ES = %.4f\n\tw2 = [%.4f, %.4f],  ES = %.4f\n", ...
    w1_1(1),w1_1(2),ES(1,1),w1_2(1),w1_2(2),ES(1,2))

fprintf("\nPtf2 VaR = 28 p = %.2f:\n", p2)
fprintf("\tw1 = [%.4f, %.4f], ES = %.4f\n\tw2 = [%.4f, %.4f],  ES = %.4f\n", ...
    w2_1(1),w2_1(2),ES(2,1),w2_2(1),w2_2(2),ES(2,2))

fprintf("\nPtf3 VaR = 28 p = %.2f:\n", p1)
fprintf("\tw1 = [%.4f, %.4f], ES = %.4f\n\tw2 = [%.4f, %.4f], ES = %.4f\n", ...
    w3_1(1),w3_1(2),ES(3,1),w3_2(1),w3_2(2),ES(3,2))

fprintf("\nPtf4 VaR = 28 p = %.2f:\n", p2)
fprintf("\tw1 = [%.4f, %.4f], ES = %.4f\n\tw2 = [%.4f, %.4f], ES = %.4f\n", ...
    w4_1(1),w4_1(2),ES(4,1),w4_2(1),w4_2(2),ES(4,2))

