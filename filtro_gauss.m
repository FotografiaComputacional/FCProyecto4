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
  
  