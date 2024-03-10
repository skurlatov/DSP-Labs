function x = quant(x, step_q)
%Функция для исследования шума квантования
x_q = round(x*step_q)/step_q;
e = x_q - x; %шум квантования

figure;
subplot(2,2,1);
plot(e(1:200));
title("Шум квантования");

subplot(2,2,2);
hist(e, 100);
title("Гистограмма шума квантования");

[AKF, dk] = xcorr(e, 100, "unbiased"); %АКФ шума квантования

subplot(2,2,3);
plot(dk, AKF);
title("АКФ шума квантования");

subplot(2,2,4);
pwelch(e, 256);

end