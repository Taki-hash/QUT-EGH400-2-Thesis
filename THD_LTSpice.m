
% INPUT INJECTION
V = [4.616e-4 2.392e-4 1.646e-4 1.266e-4 1.033e-4 8.753e-5 7.652e-5 6.909e-5 6.413e-5]; 

V1 = V(1);              % Fundamental
Vn = V(2:5);          % All other harmonics

% Calculate THD as a ratio and percentage
THD_ratio = norm(Vn) / V1;
THD_per = 100 * THD_ratio;

% Display result
fprintf('THD INPUT = %.3f%%\n', THD_per);


%% OUTPUT INJECTION
V = [5.104e-2 2.723e-3 1.589e-3 1.138e-3 8.756e-4 7.375e-4 6.226e-4 5.464e-4 4.936e-4]; 

V1 = V(1);              % Fundamental
Vn = V(2:5);          % All other harmonics

% Calculate THD as a ratio and percentage
THD_ratio = norm(Vn) / V1;
THD_per = 100 * THD_ratio;

% Display result
fprintf('THD OUTPUT = %.3f%%\n', THD_per);


