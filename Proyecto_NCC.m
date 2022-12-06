clear all
close all
clc

%LEEMOS LAS IMAGENES 
img = imread('Image_3.bmp');
temp = imread('Template.bmp');
%LAS CAMBIAMOS A ESCALAS DE GRISES
img_g = rgb2gray(img);
temp_g = rgb2gray(temp);

%OBTENEMOS SU ALTURA Y SU ANCHO
[img_H,img_W] = size(img_g);
[temp_H,temp_W] = size(temp_g);

%Limites 
xl = [1; 1];
xu = [img_W - temp_W; img_H - temp_H];

G = 250; %GENERACIONES / ITERACIONES 
N = 100; %INDIVIDUOS
D = 2; %DIMENSION 

F = 0.6; %1.2; %FACTOR DE AMPLIFICACION 
CR = 0.9; %0.6; %CONSTANTE DE RECOMBINACION 

% F = 1.2;
% CR = 0.6;

x = zeros(D,N); %SOLUCIONES
fitness = zeros(1,N); %EVALUACION DE LA FUNCION OBJETIVO
%SIEMPRE QUE TENGAMOS UN CAMBIO EN LA POSICION, ACTUALIZAMOS FITNESS	
fx_plot = zeros(1,G);

%INICIALIZACION
for i=1:N
	x(:,i) = xl+(xu-xl).*rand(D,1); %Inicializacion
 	fitness(i) = NCC(img_g,temp_g, floor(x(1,i)),floor(x(2,i)));
end

tic

for n=1:G
	%Plot_Contour(f,x,xl,xu) %Grafica2D
	for i=1:N
        %MUTACION
		I = randperm(N);
        I(I == i) = [];
        r1 = I(1);
        r2 = I(2);
        r3 = I(3);
        r4 = I(4);
        r5 = I(5);

		%DE/rand/1/bin
 		v = x(:,r1) + F*(x(:,r2) - x(:,r3));
		
		%DE/best/1/bin
 		%[~, best] = min(fitness);
 		%v = x(:,best) + F*(x(:,r2) - x(:,r3));

 		%DE/rand/2/bin
		%v = x(:,r1) + F*(x(:,r2) - x(:,r3)) + F*(x(:,r4) - x(:,r5));
		
		%DE/current to rand/1/bin
		%v = x(:,i) + F*(x(:,r1) - x(:,i)) + F*(x(:,r2) - x(:,r3));
		
		%RECOMBINACION BIN
		u = zeros(D,1);
		k = randi([1,D]);

		for j=1:D
			if rand() <= CR || j==k
				u(j) = v(j);
			else
				u(j) = x(j,i);
			end			
 		end

		
		%RECOMBINACION EXP
% 		u = x(:,i); %VECTOR DE PRUEBA
% 		j = randi([1,D]);
% 		L = 1;
% 
% 		while rand() <= CR && L <= D
% 			u(j) = v(j);
% 			j = 1 + mod(j,D);
% 			L = L+1;	
% 		end
		
		%SELECCION

        %METODO 3 / PENALIZACION
		for j=1:D
			if xl(j)<u(j) && u(j)<xu(j)
                %NADA AQUI
            else
                u(j) = xl(j) +(xu(j)-xl(j))*rand();
            end
        end

		fu = NCC(img_g,temp_g, floor(u(1)),floor(u(2))); %METODO 3

		if fu > fitness(i)
			x(:,i) = u; 
			fitness(i) = fu;
        end
    end

	[fx_plot(n), ~] = max(fitness); 
end

toc

[~, igb] = max(fitness);
%MOSTRAMOS VALORES
xp = x(1,igb);
yp = x(2,igb);

disp(['Posición en x: ' num2str(xp)]);
disp(['Posicion en y: ' num2str(yp)]);
disp(['NCC: ' num2str(fitness(igb))]);

%GRAFICA DE CONVERGENCIA
figure
hold on
grid on
plot(fx_plot, 'b-', 'LineWidth', 2)
title('Grafica de Convergencia');
xlabel('iteracion')
ylabel('f(x)');

%IMAGEN
figure
hold on
imshow(img)

line([xp xp+temp_W], [yp yp],'Color','g','LineWidth',3);
line([xp xp], [yp yp+temp_H],'Color','g','LineWidth',3);
line([xp+temp_W xp+temp_W], [yp yp+temp_H],'Color','g','LineWidth',3);
line([xp xp+temp_W], [yp+temp_H yp+temp_H],'Color','g','LineWidth',3);
