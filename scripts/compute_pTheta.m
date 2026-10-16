
function pTheta = compute_pTheta(T, dx)
  % Compute the p-Theta edge detection attribute.

  [dTdx, dTdz] = gradient(T, dx);

  epsilon = 1e-6;

  pTheta = atan(dTdz ./ (dTdx + epsilon));
end
