function res=filtro_gausss(im,S)

  L = round(S*2); 
  Gs = fspecial('gaussian',2*L+1,S); %Creación máscara gaussiana
  im=double(im); 
  
  [N,M,P]=size(im);
  
  % Añadimos margen de L filas y columnas por cada lado
  % repitiendo la 1ª/última fila y la 1ª/última columna
  % Equivale a opción symmetric en imfilter
  im=[im(:,L:-1:1,:) im im(:,M:-1:M-L+1,:)];
  im=[im(L:-1:1,:,:);im;im(N:-1:N-L+1,:,:)];
      
  res=im*0; s=(-L:L); 
  for k=L+1:L+N
    for j=L+1:L+M               
        vec = im(k+s,j+s,:);   
        vec = vec.*Gs;    
        res(k,j,:) = sum(sum(vec));
    end
  end
  
  % Quitamos el margen que le hemos añadido al principio
  res = res(L+(1:N),L+(1:M),:);
  end
  
function show_detail(det)
    det_abs = abs(det);
    I = 0.3 * det_abs(:,:,1) + 0.55 * det_abs(:,:,2) + 0.15 * det_abs(:,:,3);
    imagesc(I);
    colormap('hot');
    colorbar('vert');
    axis image;
    title('Detalle visualizado');
  end

im = imread('img1.jpg');
im = double(im);

S = 2;
imS = filtro_gausss(im, S);
detalle = im - imS;
show_detail(detalle);

% Realce del detalle
alpha = 1.25; 
umbral = 3.0;

figure;
mascara = abs(detalle) > umbral;
im_realzada_umbral = im + alpha * (detalle .* mascara);
imshow(uint8(im_realzada_umbral));
title('Realce del detalle con umbral');

% Imagen mitad original + mitad realzada
figure;
compuesta = im;
im_realzada = im + alpha * detalle;
compuesta(701:end,:,:) = im_realzada(701:end,:,:);
compuesta(700,:,:) = 0;  
imshow(uint8(compuesta));
title('Imagen compuesta: mitad original + mitad realzada');

% ------------------------------------------------------------
% Parte 2
function res = filtro_bilat(im, S, R)
  L = round(2*S);  
  Gs = fspecial('gaussian', 2*L+1, S);  
  im = double(im); 

  [N, M, P] = size(im);
  im = [im(:,L:-1:1,:) im im(:,M:-1:M-L+1,:)];
  im = [im(L:-1:1,:,:); im; im(N:-1:N-L+1,:,:)];
  
  res = zeros(size(im));  % inicializamos resultado
  s = -L:L;

  for k = L+1 : L+N
    for j = L+1 : L+M
      vec = im(k+s,j+s,:);  %Extraemos la vecindad
      centro = im(k,j,:);   %Pixel central
      D = vec - centro;     %Diferencias entre el pixel central y los de la vecindad
      D = D / R;            %Division por el ancho
      D = sum(D.^2, 3);    %Sumamos los valores de D^2 en su tercera dimension
      
      Gr= exp(-0.5 *D);    %Obtener los coeficientes de la 2a gaussiana
      G = Gr .* Gs;         %Multiplicamos punto a punto la nueva mascara con la anterior
      G = G / sum(G(:));    

      for i = 1 : P
        patch = vec(:,:,i);     %Extraemos la vecindad del canal i
        res(k,j,i) = sum(sum(patch .* G)); %Aplicamos la mascara
      end
    end
  end

  res = res(L+(1:N),L+(1:M),:);
end

im = imread('img1.jpg');
im = double(im);
im = imresize(im, 1/4);

R = 20;
L = round(2*S);  
bf1 = filtro_bilat(im, S, R);
bf2 = imbilatfilt(im,R^2,S,'Padding','symmetric','NeighborhoodSize',2*L+1); 

diff = bf1 - bf2;
max_diff = max(diff(:));
min_diff = min(diff(:));

disp(['Máximo error: ', num2str(max_diff)]);
disp(['Mínimo error: ', num2str(min_diff)]);

%Aplicamos el filtro bilateral a la imagen original
im = imread('img1.jpg');
im = double(im);

im_bilart = filtro_bilat(im, S, R);
detalle_bilat = im - im_bilart;
figure;
show_detail(detalle_bilat);
title('Detalle de la imagen original con filtro bilateral');

% -----------------------------
im = imread('img1.jpg');
im = double(im);
% 1. Suavizado con filtro gaussiano 10 veces
S = 4;

gauss_result = im;
for i = 1:10
    gauss_result = imgaussfilt(gauss_result, S, 'Padding', 'symmetric');
end

figure;
imshow(uint8(gauss_result));
title('10 filtrados gaussianos (S = 4)');

% -----------------------------
% 2. Suavizado con filtro bilateral 10 veces
R = 20;
L = round(2*S);

bilat_result = im;
for i = 1:10
    bilat_result = imbilatfilt(bilat_result, R^2, S, 'Padding', 'symmetric');
end

figure;
imshow(uint8(bilat_result));
title('10 filtrados bilaterales (S = 4, R = 20)');

% -----------------------------
% Cross bilateral filter
function res=filtro_cross(im, im2, S, R)
  L = round(2*S);  
  Gs = fspecial('gaussian', 2*L+1, S);  
  im = double(im); 
  im2 = double(im2);

  [N, M, P] = size(im);
  im = [im(:,L:-1:1,:) im im(:,M:-1:M-L+1,:)];
  im = [im(L:-1:1,:,:); im; im(N:-1:N-L+1,:,:)];
  im2 = [im2(:,L:-1:1,:) im2 im2(:,M:-1:M-L+1,:)];
  im2 = [im2(L:-1:1,:,:); im2; im2(N:-1:N-L+1,:,:)];
  
  res = zeros(size(im));  % inicializamos resultado
  s = -L:L;

  for k = L+1 : L+N
    for j = L+1 : L+M
      vec = im(k+s,j+s,:);  %Extraemos la vecindad
      vec2 = im2(k+s,j+s,:); %Extraemos la vecindad de la imagen 2
      centro = im2(k,j,:);   %Pixel central

      D = vec2 - centro;     %Diferencias entre el pixel central y los de la vecindad
      D = D / R;            %Division por el ancho
      D = sum(D.^2, 3);    %Sumamos los valores de D^2 en su tercera dimension
      
      Gr= exp(-0.5 *D);    %Obtener los coeficientes de la 2a gaussiana
      G = Gr .* Gs;         %Multiplicamos punto a punto la nueva mascara con la anterior
      G = G / sum(G(:));

      for i = 1 : P
        patch = vec(:,:,i);     %Extraemos la vecindad del canal i
        res(k,j,i) = sum(sum(patch .* G)); %Aplicamos la mascara
      end
    end
  end

  res = res(L+(1:N),L+(1:M),:);
end

im_flash = imread('flash.jpg');
im_noflash = imread('no_flash.jpg');

im_flash = double(im_flash);
im_noflash = double(im_noflash);

% Filtro gaussiano a la imagen sin flash
S = 5;
noflash_gauss = imgaussfilt(im_noflash, S, 'Padding', 'symmetric');
figure;
imshow(uint8(noflash_gauss));
title('No flash suavizada con filtro gaussiano');

%Aplicamos el filtro cross bilateral
R = 10;  
noflash_cross = filtro_cross(im_noflash, im_flash, S, R);
figure;
imshow(uint8(noflash_cross));
title('No flash suavizada con filtro cross-bilateral');

%Juntamos para comparar
comparacion = [uint8(im_noflash), uint8(noflash_gauss), uint8(noflash_cross)];
figure;
imshow(comparacion);
title('Original | Gaussiano | Cross-bilateral');

