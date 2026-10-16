
function ISA = compute_ISA(T, dx)
  % Compute the Inclination of the Source Analytic Signal (ISA).

  [dTdx, dTdz] = gradient(T, dx);

  ASA = sqrt(dTdx.^2 + dTdz.^2);
  epsilon = 1e-6;

  ISA = atan(dTdz ./ (ASA + epsilon));
end
