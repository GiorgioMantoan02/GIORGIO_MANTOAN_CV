clear all
clc 
close all

%% ES 1 Apple

filename = 'Gruppo 14 Lab 3.xlsx';
sheetname= 'Esercizi aula 1';
range = 'A:B';
data = readmatrix(filename, 'Range', range, 'Sheet', sheetname);

% tengo solo le date che mi servono
data = data(4:end,2);

returns = log(data(1:end-1)./data(2:end)); % vicino --> lontano
returns = flipud(returns); % riordinamento

figure(1)

plot(returns,'r')
ylim([-0.15,0.15])
xlabel('days')
ylabel('log-returns')
title('plot log-returns Apple')

%% ES 2 istogramma

% istogramma
figure(2)
hist(returns)

% media
mu = mean(returns);

% deviazione standard 
sigma = std(returns);

% Test di Kolmogorov-Smirov rispetto a una distribuzione normale standard
% Se p-value <= 0.05, rifiuta ipotesi nulla
% --> i dati non sono normalmente distribuiti

[h_ks, p_ks] = kstest((data - mu)/sigma);
disp(['Kolmogorov-Smirov p-value: ', num2str(p_ks)]);

figure(3)
qqplot((data-mu)/sigma)

%% Curtosi

skrew = skewness (returns);
curtosi = kurtosis(returns);
excess_kurt = curtosi - 3;

disp(['eccesso di curtosi: ', num2str(excess_kurt)])


%% FRONTIERA DEI PORTAFOGLI

filename = 'Gruppo 14 Lab 3.xlsx';
sheetname= 'Esercizi aula 2';
range = 'A:D';
data_FP = readmatrix(filename, 'Range', range, 'Sheet', sheetname);

%APPLE

data_apple = data_FP(5:end,2);  
returns_apple = log(data_apple(1:end-1)./data_apple(2:end));
returns_apple = flipud(returns_apple);

mu_apple = mean(returns_apple);

%AMAZON

data_amazon = data_FP(5:end,3);
returns_amazon = log(data_amazon(1:end-1)./data_amazon(2:end));
returns_amazon = flipud(returns_amazon);

mu_amazon = mean(returns_amazon);

%FACEBOOK

data_facebook = data_FP(5:end,4);
returns_facebook = log(data_facebook(1:end-1)./data_facebook(2:end));
returns_facebook = flipud(returns_facebook);

mu_facebook = mean(returns_facebook);


%calcolo di w

M = [returns_apple, returns_amazon, returns_facebook];  % matrice dei rendimenti

V = cov(M);                                             % matrice delle covarianze
e = [mu_apple, mu_amazon, mu_facebook]';                % vettore delle medie

One = ones(size(e)); 
invV= inv(V);

A = One'*invV*e;
B = e'*invV*e;
C = One'*invV*One;
D = B*C-A^2;

g = (B*(invV*One)-A*(invV*e))/D;

h = (-A*(invV*One)+C*(invV*e))/D;


w = @(R) g+h*R;     % w in funzione di E(r):=R

Varianza = @(R) (C*R.^2-2*A*R+B)/D;

R = linspace(min(e)-0.002, max(e)+0.002, 1000);

figure(4)

plot(Varianza(R), R,'LineWidth', 2,'Color','b');
xlabel('Varianza'); ylabel('Return');
hold on
grid on
plot(diag(V), e,'*', 'LineWidth', 1, 'Color', 'r')
legend ('Frontiera Portafogli','Amazon - Facebook - Apple')

hold off


%% (1) Portafoglio A: AMAZON + FACEBOOK

% Rendimenti, varianza e media

M_A = [returns_amazon, returns_facebook];
V_A = cov(M_A);
e_A = [mu_amazon, mu_facebook]';

% Frontiera portafoglio efficiente 

%Chiamo "A" -> "1" per comodità nei parametri

One1 = ones(size(e_A));
invV1 = inv(V_A);

A1 = One1' * invV1 * e_A;
B1 = e_A' * invV1 * e_A;
C1 = One1' * invV1 * One1;
D1 = B1*C1 - A1^2;

g1 = (B1*(invV1*One1) - A1*(invV1*e_A)) / D1;
h1 = (-A1*(invV1*One1) + C1*(invV1*e_A)) / D1;

% Pesi e varianza in funzione del rendimento E(r_A):=R1
w1 = @(R1) g1 + h1 * R1;
Varianza1 = @(R1) (C1*R1.^2 - 2*A1*R1 + B1) / D1;

% Valori di R1 utili per il grafico
R1 = linspace(min(e_A)-0.001, max(e_A)+0.001, 1000);

figure (5);
plot(Varianza1(R1), R1, 'b', 'LineWidth', 2);
title('Frontiera A: Amazon - Facebook');
xlabel('Varianza');
ylabel('Rendimento Atteso');
grid on;
hold on;
plot(diag(V_A), e_A, '*', 'LineWidth', 2, 'color', 'r');
legend('FP', 'Amazon - Facebook');


%% (2) portafoglio B: APPLE + FACEBOOK

% Rendimenti, varianza e media

M_B = [returns_apple, returns_facebook];
V_B = cov(M_B);
e_B = [mu_apple; mu_facebook];

% Frontiera portafoglio efficiente

%Chiamo "B" -> "2" per comodità nei parametri

One2 = ones(size(e_B));
invV2 = inv(V_B);

A2 = One2' * invV2 * e_B;
B2 = e_B' * invV2 * e_B;
C2 = One2' * invV2 * One2;
D2 = B2*C2 - A2^2;

g2 = (B2*(invV2*One2) - A2*(invV2*e_B)) / D2;
h2 = (-A2*(invV2*One2) + C2*(invV2*e_B)) / D2;

% Pesi e varianza in funzione del rendimento E(r_B):=R2
w2 = @(R2) g2 + h2 * R2;
Varianza2 = @(R2) (C2*R2.^2 - 2*A2*R2 + B2) / D2;

% Valori di R2 utili per il grafico
R2 = linspace(min(e_B)-0.001, max(e_B)+0.001, 1000);

figure (6);
plot(Varianza2(R2), R2, 'b', 'LineWidth', 2);
title('Frontiera B: Apple - Facebook');
xlabel('Varianza');
ylabel('Rendimento Atteso');
grid on;
hold on;
plot(diag(V_B), e_B, '*', 'LineWidth', 2,'Color', 'r');
legend('FP', 'Apple - Facebook');


%% (3) Portafoglio MVP per A, B, tutti i titoli

w_minVar_A = invV1 * One1 / C1;

fprintf('\nPortafoglio a varianza minima (Amazon - Facebook):\n');
fprintf('w_Amazon:   %.4f\n', w_minVar_A(1));
fprintf('w_Facebook: %.4f\n', w_minVar_A(2));

w_minVar_B = invV2 * One2 / C2;

fprintf('\nPortafoglio a varianza minima (Apple - Facebook):\n');
fprintf('w_Apple:    %.4f\n', w_minVar_B(1));
fprintf('w_Facebook: %.4f\n', w_minVar_B(2));

w_minVar_all = invV * One / C;

fprintf('\nPortafoglio a varianza minima (Amazon - Facebook):\n');
fprintf('w_Apple:    %.4f\n', w_minVar_all(1));
fprintf('w_Amazon:   %.4f\n', w_minVar_all(2));
fprintf('w_Facebook: %.4f\n', w_minVar_all(3));


%% (4) Portafoglio ottimo con utilità quadratica

% Parametri
rf = 1.01;
rf_log = log(rf);
b = 0.001;
m = 252*e;
V_anno = 252*V;
invV_anno = inv(V_anno);

exc_ret = m - rf_log * One; % excess return

w_opt = (((1-b*rf)/b)*invV_anno*exc_ret)/(1+exc_ret'*invV_anno*exc_ret);


% Calcolo rendimento e varianza portafoglio
R_p = rf + exc_ret' * invV_anno * exc_ret / b;
Var_p = (1/b^2) * (exc_ret' * invV_anno * exc_ret);

fprintf('\nPortafoglio ottimo:\n');
fprintf('w_Apple:    %.4f\n', w_opt(1));
fprintf('w_Amazon:   %.4f\n', w_opt(2));
fprintf('w_Facebook: %.4f\n', w_opt(3));
fprintf('Rendimento atteso: %.4f\n', R_p);
fprintf('Varianza: %.4f\n', Var_p);

%% (5) Portafoglio ottimo con utilità esponenziale
% Parametri
a = 0.15;
rf = 1.01;
x = 1;

sigma2_apple = V_anno(1,1);
sigma2_amazon = V_anno(2,2);
sigma2_facebook = V_anno(3,3);

w_apple= (m(1)-rf_log)/(a*sigma2_apple);
w_amazon= (m(2)-rf_log)/(a*sigma2_amazon);
w_facebook= (m(3)-rf_log)/(a*sigma2_facebook);

E_apple = (-1/a)*exp(-a*(w_apple*m(1)+(x-w_apple)*rf_log)+a^2*w_apple^2*sigma2_apple/2);
E_amazon = (-1/a)*exp(-a*(w_amazon*m(2)+(x-w_amazon)*rf_log)+a^2*w_amazon^2*sigma2_amazon/2);
E_facebook = (-1/a)*exp(-a*(w_facebook*m(3)+(x-w_facebook)*rf_log)+a^2*w_facebook^2*sigma2_facebook/2);

% Titolo migliore
[~, idx] = max([E_apple, E_amazon, E_facebook]);
titoli = {'Apple', 'Amazon', 'Facebook'};
disp(['Miglior titolo da affiancare al risk-free: ', titoli{idx}]);