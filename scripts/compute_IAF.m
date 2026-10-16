function IAF = compute_IAF(T, dx)
  % Compute the Inclination-Based Angular Fusion (IAF) attribute.
  % Manuscript definition: IAF = (ISA + Theta) / 2.
  % ISA = atan(Tz / (ASA + epsilon))
  % Theta = acos(Tz / (ASA + epsilon))

  [dTdx, dTdz] = gradient(T, dx);
  ASA = sqrt(dTdx.^2 + dTdz.^2);
  epsilon = 1e-6;

  ratio = dTdz ./ (ASA + epsilon);
  ratio = max(-1, min(1, ratio));  % numerical safety

  ISA = atan(ratio);
  Theta = acos(ratio);
  IAF = (ISA + Theta) / 2;
end
