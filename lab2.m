% Лабораторная работа №2 -  ДИСКРЕТНЫЕ ФИЛЬТРЫ
% Данная функция реализует все пункты работы
function filter_analysis = lab2(U1, U2, U3, U4, T1, T2, Fd,a1,b1)

% Создание дискретного сигнала из первой лабораторной
T = 0:(1/Fd):T2; %Общая ось времени
signal = zeros([size(T)]);
for i = 1:length(T)
    if T(i) < T1
       signal(i) = T(i)*(U2 - U1)/T1 + U1; %из канонического уравнения прямой
    else
       signal(i) = (T(i)-T1)*(U4 - U3)/(T2-T1) + U3;
    end
end
signal = [signal zeros(size(signal))];

% Создание вектора для максимальных
% по модулю значений сигналов при разных формах реализации
modmax = zeros(3);

% Дискретный фильтр (прямая форма)
sig_out1 = filter(b1,a1,signal);
modmax(1) = max(abs(sig_out1));
% Дискретный фильтр (каноническая форма)
sig_out2 = filter(1,a1,signal); % Получение значений сигнала в элементах памяти
modmax(2) = max(abs(sig_out2));
% Дискретный фильтр (транспонированная форма)
states = []; % Матрица внутренних состояний фильтра
cur = []; % Текущее состояние фильтра
for k = 1:length(signal)
    [sig_out3(k),cur] = filter(b1,a1,signal(k),cur);
    states = [states cur];
end
modmax(3) = max(max(abs(sig_out3)));
% Получение аналитического выражения ИХ фильтра
% Нахождение векторов вычетов, полюсов и целой части функции передачи
[r,p,k] = residuez(b1,a1);
% Модули и фазы полюсов и вычетов
r_abs = abs(r);
p_abs = abs(p);
r_ang = angle(r);
p_ang = angle(p);

%Построение всех графиков лабораторной работы

figure;
% Входной сигнал из первой лабораторной
subplot(1,3,1);
stem(signal);
% Выходной сигнал дискретного фильтра прямой формы
subplot(1,3,2);
stem(sig_out1);
% Сигнал дискретного фильтра канонической формы в элементах памяти
subplot(1,3,3);
stem(sig_out2);
% Графики внутренних состояний дискретного фильтра транспонированной формы
figure;
plot(states');
% Графики характеристик фильтра
% Все нужные изображения копируются вручную!
% Отсюда берутся:
% АЧХ, ФЧХ, групповая задержка, ИХ, график нулей и полюсов
fvtool(b1,a1);
end