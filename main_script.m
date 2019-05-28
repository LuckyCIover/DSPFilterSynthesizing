%% Данные варианта

M = 5;

pass_band = [0 0.6]*pi;
supp_band = [0.8 1]*pi;

%% задание 1

j = 0:1:M;
w = pi * (j + 0.5) / (M + 1);

pass_band_indices  = find(w >= pass_band(1) & w <= pass_band(2));
trans_band_indices = find(w >  pass_band(2) & w <  supp_band(1));
supp_band_indices  = find(w >= supp_band(1) & w <= supp_band(2));

Kd_p =  ones(1, length(pass_band_indices));
Kd_t = zeros(1, length(trans_band_indices));
Kd_s = zeros(1, length(supp_band_indices));

Kd = [Kd_p Kd_t Kd_s];

b = find_b(Kd, w, M);
[AFC, AFC1] = build_AFC(b, M);

% ----- Линейная шкала ----- %
figure, hold on, grid on

% AЧХ
plot(w, Kd, 'b*');
fplot(AFC, [-pi, 2*pi], 'b', 'LineWidth', 1);

% идеальная АЧХ
IAFC = @(w) 1 - heaviside(w - (pass_band(2) + supp_band(1))/2);
fplot(IAFC, [-pi, 2*pi], 'm', 'LineWidth', 1);

title('АЧХ')
xlabel('w'), ylabel('|K(w)|')
legend('w_j', 'АЧХ', 'Идеальная АЧХ')

% ----- Логарифмическая шкала ----- %
figure, hold on, grid on

w_ = -pi:pi/50:2*pi;
AFC_ = AFC(w_);
AFC_(AFC_ < 10^-7) = 10^-7;
AFC_log = 20*log10(AFC_);

% AЧХ
plot(w_, AFC_log, 'b', 'LineWidth', 1);

title('АЧХ')
xlabel('w'), ylabel('20lg|K(w)/C|')

% ----- ФЧХ ----- %
figure, hold on, grid on
FFC = angle(exp(-i*w_*M).*AFC1(w_));
plot(w_, FFC, 'b', 'LineWidth', 1);
title('ФЧХ')
xlabel('w'), ylabel('фаза')

%% Задание 2

synthesize_filter([0 0.6]*pi, [0.8 1]*pi, 5, 1);

%% Задание 3

M_ = 9;

pass_err = 1;
supp_err = 1;
while pass_err > 0.0125 && supp_err > 0.015
    M_ = M_ + 1;
    [b, ~, pass_err, supp_err] = synthesize_filter([0 0.6]*pi, ...
        [0.8 1]*pi, M_);
end

fprintf('Удовлетворяет требованиям при M = %i', M_);

%% Задание 4

B = [b flip(b(1:end-1))];

n = 1:50;
X1 = sin(0.3*pi*n);

Y1 = filter(B, 1, X1);

figure, hold on, grid on

plot(n, X1, 'b--', 'LineWidth', 1)
plot(n, Y1, 'r', 'LineWidth', 1)

%% Задание 5

B = [b flip(b(1:end-1))];

img = double(imread('var9.png'));

figure
imshow(uint8(img))

for i = 1:1:size(img,1)
    img(i,:) = filter(B, 1, img(i,:));
end

img(:,1:end-M) = img(:,M+1:end);

for i = 1:1:size(img,2)
    img(:,i) = filter(B, 1, img(:,i));
end

img(1:end-M,:) = img(M+1:end,:);

figure
imshow(uint8(img))

%% вспомогательные функции

function[b] = find_b(Kd, w, M)
    b = zeros(1, length(w));

    for k = 0:1:M
        S = 0;
        for j = 0:1:M
            S = S + Kd(j+1) * cos(k*w(j+1));
        end
        b(M-k+1) = S / (M + 1);
    end
end

function [abs_AFC, AFC] = build_AFC(b, M)

    syms w
    
    S = 0;
    for k = 1:1:M 
        S = S + b(M-k+1)*cos(w*k);
    end

    abs_AFC = matlabFunction(abs(b(M+1) + 2*S)); 
    AFC = matlabFunction(b(M+1) + 2*S);
    
end
