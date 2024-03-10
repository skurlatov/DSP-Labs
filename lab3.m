%%ТРЕТЬЯ ЛАБОРАТОРНАЯ.
function x = lab3(U1, U2, U3, U4, T1, T2, Fd)
    %Формирование сигнала
    T = 0;
    Sd = 0;
    [T, Sd] = signal_form(U1, U2, U3, U4, T1, T2, Fd);
    spectrum_fft = fft(Sd); %расчет спектральных отсчётов
    E0 = sum(abs(Sd).^2); %расчёт энергии исходного сигнала
    %Оценка ширины спектра
    N_max = 0;
    E0_n = 0;
    spectrum_n = spectrum_fft;
    while (E0_n < 0.9*E0)
        spectrum_n = spectrum_fft;
        spectrum_n(N_max+1:length(spectrum_n)-N_max+1) = 0;
        E0_n = sum(abs(ifft(spectrum_n)).^2);
        N_max = N_max+1;
    end

    %Дополнение нулями
    S_zeros = [Sd zeros([size(Sd)])];
    spectrum_zeros = fft(S_zeros);

    %Измерение скорости расчётов по теории и через fft
    N = [64 128 256 512 1024 2048 4096 8192]; % вектор размеров ДПФ
    K = 750; %количество циклов расчета. Изменяемый параметр
    T_N_theory = zeros([size(N)]); %затраченное время на расчет по прямой формуле
    T_N_fft = zeros([size(N)]); %затраченное время на расчет по быстрому алгоритму
    for i = 1:length(N)
        % дополнение сигнала нулями до длины N
        x1 = [Sd zeros(1, N(i)-length(Sd))];
        D = dftmtx(N(i)); % матрица ДПФ
        y = zeros(1, N(i)); % массив для результатов ДПФ
        tic % старт таймера
        for k = 1:K % цикл для измерения времени
            x1 * D; % вычисление ДПФ по прямой формуле
        end
        T_N_theory(i) = toc; %запись во временной массив
        tic
        for k = 1:(K*300) %K*300 - изменяемый параметр
            fft(Sd, N(i));
        end
        T_N_fft(i) = toc;
    end
    T_one_theory = T_N_theory ./ K; %время на однократное вычисление ДПФ по прямой формуле
    T_one_fft = T_N_fft ./ (K*300); %время на однократное вычисление ДПФ по быстрому алгоритму

    %Графики и выводы
    E0
    N_max
    K
    K*300
    k1 = 0.9
    k2 = 0.4
    N
    T_N_theory
    T_one_theory*10^9
    T_N_fft
    T_one_fft*10^9
    
    figure;

    subplot(4, 2, 1);
    plot(T, Sd);
    title("График дискретного сигнала (plot)"); xlabel("t, мс"); ylabel("s, В");

    subplot(4, 2, 2);
    stem(abs(spectrum_fft));
    title("Модуль FFT"); xlabel("N, отсчеты"); ylabel("s, В/Гц");

    subplot(4, 2, 3);
    stem(angle(spectrum_fft));
    title("Фазы FFT"); xlabel("N, отсчеты"); ylabel("\phi, рад");

    subplot(4,2,4);
    stem(T, Sd);
    hold on;
    stem(T, ifft(spectrum_n));
    title("Исходный сигнал и сигнал с 0.9 E0"); xlabel("T, отсчеты"); ylabel("s, В");

    subplot(4, 2, 5);
    stem(abs(spectrum_zeros));
    title("Модуль FFT сигнала, дополненного нулями"); xlabel("N, отсчеты"); ylabel("s, В/Гц");
    
    subplot(4, 2, 6);
    stem(angle(spectrum_zeros));
    title("Фазы FFT сигнала, дополненного нулями"); xlabel("N, отсчеты"); ylabel("\phi, рад");

    subplot(4, 2, 7);
    plot(log2(N), log(T_one_theory), log2(N), log((k1 * N.^2)*10^-9))
    title("Время, затраченное на однократное вычисление ДПФ прямым способом"); xlabel("{log}_2 (N)"); ylabel("log(t)");

    subplot(4, 2, 8);
    plot(log2(N), log(T_one_fft), log2(N), log((k2 * N .* log2(N))*10^-9))
    title("Время, затраченное на однократное вычисление ДПФ быстрым способом"); xlabel("{log}_2 (N)"); ylabel("log(t)");
end
