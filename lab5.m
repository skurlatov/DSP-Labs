%ПЯТАЯ ЛАБОРАТОРНАЯ
function x = lab5(A, w0, phi0, Type, fcp, Apass, Astop)
load mtlb;
k = 1:(10^5-1);
signal = A*cos(w0*k+phi0);
WGN = randn(1, 10^5);
quant(signal, 256); %для исходного сигнала
quant(signal, 16); %для того же сигнала, но с шагом квантования 1/16
quant(signal, 1); %для того же сигнала, но с шагом квантования 1
quant(WGN/max(abs(WGN)), 256); %для БГШ
quant(mtlb/max(abs(mtlb)), 256); %для речевого сигнала

fdatool; %3-5 пункты выполняются тут. В файлах приложены полученные графики и файл фильтра для нашего варианта

figure;
w=0:pi/100:pi;
A=5*0.4*10^(-5);%масштабный коэффициент
%теоретическая СПМ, коэффициенты взяты из результатов выполнения 3-5 пунктов
Wn=A./(abs(1-3.572*exp(-1i*w)+4.808*exp(-1i*w*2)-2.886*exp(-1i*w*3)+0.0521*exp(-1i*w*4))).^2;
plot(w,10*log10(Wn));
end