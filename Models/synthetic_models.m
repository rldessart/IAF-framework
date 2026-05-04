function [T, dx] = synthetic_models(model_type)

  % ============================================
  % GRID
  % ============================================
  nx = 100;
  nz = 100;
  dx = 1;

  [X,Z] = meshgrid(1:nx,1:nz);

  % ============================================
  % INITIALIZE MODEL
  % ============================================
  T = zeros(nz,nx);

  % ============================================
  % SELECT MODEL
  % ============================================
  switch model_type

    case 1  % circular body
      T((X-50).^2 + (Z-50).^2 <= 15^2) = 1;

    case 2  % MULTI-SOURCE + NOISE

  % Dique principal
  T(:,48:52) = 1;

  % Corpo circular 1
  T((X-30).^2 + (Z-70).^2 <= 12^2) = 1.2;

  % Corpo circular 2
  T((X-70).^2 + (Z-40).^2 <= 10^2) = 0.8;

  % Corpo adicional
  T((X-60).^2 + (Z-80).^2 <= 6^2) = 1;

  % Suavização
  G = gaussian_kernel(15,4);
  T = conv2(T, G, 'same');

  % Ruído
  T = T + 0.05 * randn(size(T));

    case 3  % vertical contact
      T(:,50:end) = 1;

% leve suavização adicional
T = conv2(T, gaussian_kernel(21,5), 'same');

    otherwise
      error('Invalid model type')
  end

  % ============================================
  % SMOOTHING (simulate magnetic response)
  % ============================================
  G = gaussian_kernel(9, 2);
  T = conv2(T, G, 'same');

end


% ============================================
% GAUSSIAN KERNEL
% ============================================
function G = gaussian_kernel(sz, sigma)
  [x,y] = meshgrid(-floor(sz/2):floor(sz/2));
  G = exp(-(x.^2 + y.^2)/(2*sigma^2));
  G = G / sum(G(:));
end
