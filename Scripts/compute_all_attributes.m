function [ISA, pTheta, IAF] = compute_all_attributes(T, dx)

  % Derivatives
  [dTdx, dTdz] = gradient(T, dx);

  % Magnitude of gradient
  grad_mag = sqrt(dTdx.^2 + dTdz.^2);

  epsilon = 1e-6;

  % --------------------------------------------
  % ISA
  % --------------------------------------------
  ISA = atan(dTdz ./ (grad_mag + epsilon));

  % --------------------------------------------
  % p-Theta
  % --------------------------------------------
  pTheta = atan(dTdz ./ (dTdx + epsilon));

  % --------------------------------------------
  % IAF (melhorado)
  % Combina amplitude + direção
  % --------------------------------------------
  IAF = atan( (dTdz .* grad_mag) ./ (dTdx.^2 + dTdz.^2 + epsilon) );

end
