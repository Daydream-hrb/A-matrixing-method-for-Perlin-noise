
clear;
clc;
close;

% Config
noise_size = 256;
i_min = 2;
i_max = 6;
Grad_factor = 0.2;
Polarization_factor = 1;

% Fractal sum
perlin_noise = zeros(noise_size, noise_size, 'single');
tic
for i = i_min: i_max
    perlin_noise = perlin_noise + perlin(noise_size, 2 ^ i);
    toc
end

% Erosion (you can process the raw noise anyway you like)
% Scale normalization
noise_max = max(perlin_noise, [], 'all');
noise_min = min(perlin_noise, [], 'all');
perlin_noise = (perlin_noise - noise_min) ./ (noise_max - noise_min);
% Gradient enhencement
[Gx, Gy, Gz] = gradient(perlin_noise);
G = sqrt(Gx .^ 2 + Gy .^ 2 + Gz .^ 2);
perlin_noise = perlin_noise .* (1 ./ (1 + Grad_factor .* G));
% Polarization
perlin_noise = perlin_noise .^ 1.2;

% Colors and contour value lists
color_list = [
    213, 216, 234;
    239, 199, 81;
    227, 130, 192;
    104, 78, 168;
    130, 23, 34] ./ 255;
% contour_list is adjustable as long as length < 5 and values in [0, 1]
% contour_list = [0.1, 0.3, 0.5, 0.7, 0.9];
contour_list = [0.2, 0.4, 0.6, 0.8];

% Visualization
figure('Color', 'w')
[X, Y, Z] = ndgrid(1: noise_size, 1: noise_size, 1: noise_size);
for i = 1: 1: size(contour_list, 2)
    s = isosurface(X, Y, Z, perlin_noise, contour_list(i));
    p = patch(s);
    set(p, 'FaceColor', color_list(i, :));  
    set(p, 'EdgeColor', 'none');
    set(p, 'FaceAlpha', 0.6);
    hold on
end
grid on
axis equal
xlabel('\it{x}')
ylabel('\it{y}')
zlabel('\it{z}')
legend(cellstr(string(contour_list)))
view(-24, 24)
set(gca, 'FontName', 'Times New Roman')

% Core function
function perlin_noise = perlin(noise_size, freq)

    % Config
    xyz_step = noise_size ./ freq;
    Amp = 10 / freq;

    % Gradient matrix of grid nodes
    grad_mat = Amp * randn(noise_size / xyz_step + 1, noise_size / xyz_step + 1, ...
        noise_size / xyz_step + 1, 3, 'single'); % n x n x n x 3
    
    % Martix of all nodes
    [X, Y, Z] = ndgrid(1: noise_size, 1: noise_size, 1: noise_size);
    
    % Unit increment & Gradient matrix of all nodes
    dx = single(mod(X - 0.5, xyz_step) / xyz_step); % n x n x n
    dy = single(mod(Y - 0.5, xyz_step) / xyz_step); % n x n x n
    dz = single(mod(Z - 0.5, xyz_step) / xyz_step); % n x n x n
    Gx = single(repelem(squeeze(grad_mat(:, :, :, 1)), xyz_step, xyz_step, xyz_step)); % n x n x n
    Gy = single(repelem(squeeze(grad_mat(:, :, :, 2)), xyz_step, xyz_step, xyz_step)); % n x n x n
    Gz = single(repelem(squeeze(grad_mat(:, :, :, 3)), xyz_step, xyz_step, xyz_step)); % n x n x n
    
    % Gradient bias list
    Gb = single([1, xyz_step + 1, xyz_step * freq, xyz_step * (freq + 1)]);
    Gb0 = Gb(1): Gb(3);
    Gb1 = Gb(2): Gb(4);

    % Integral (in the order of x->y->z)
    dot000 = Gx(Gb0, Gb0, Gb0) .* dx + Gy(Gb0, Gb0, Gb0) .* dy + Gz(Gb0, Gb0, Gb0) .* dz; % n x n x n
    dot001 = Gx(Gb0, Gb0, Gb1) .* dx + Gy(Gb0, Gb0, Gb1) .* dy + Gz(Gb0, Gb0, Gb1) .* (dz - 1); % n x n x n
    dot010 = Gx(Gb0, Gb1, Gb0) .* dx + Gy(Gb0, Gb1, Gb0) .* (dy - 1) + Gz(Gb0, Gb1, Gb0) .* dz; % n x n x n
    dot011 = Gx(Gb0, Gb1, Gb1) .* dx + Gy(Gb0, Gb1, Gb1) .* (dy - 1) + Gz(Gb0, Gb1, Gb1) .* (dz - 1); % n x n x n
    dot100 = Gx(Gb1, Gb0, Gb0) .* (dx - 1) + Gy(Gb1, Gb0, Gb0) .* dy + Gz(Gb1, Gb0, Gb0) .* dz; % n x n x n
    dot101 = Gx(Gb1, Gb0, Gb1) .* (dx - 1) + Gy(Gb1, Gb0, Gb1) .* dy + Gz(Gb1, Gb0, Gb1) .* (dz - 1); % n x n x n
    dot110 = Gx(Gb1, Gb1, Gb0) .* (dx - 1) + Gy(Gb1, Gb1, Gb0) .* (dy - 1) + Gz(Gb1, Gb1, Gb0) .* dz; % n x n x n
    dot111 = Gx(Gb1, Gb1, Gb1) .* (dx - 1) + Gy(Gb1, Gb1, Gb1) .* (dy - 1) + Gz(Gb1, Gb1, Gb1) .* (dz - 1); % n x n x n
    
    % Smoothing
    tx = fade(dx);
    ty = fade(dy);
    tz = fade(dz);
    
    % Bilinear interplotation
    perlin_noise = lerp(tz, ...
        lerp(ty, lerp(tx, dot000, dot100), lerp(tx, dot010, dot110)), ...
        lerp(ty, lerp(tx, dot001, dot101), lerp(tx, dot011, dot111)));

end

% Smoothing func
function t_ = fade(t)
    t_ = t .* t .* t .* (t .* (t .* 6 - 15) + 10);
end

% Bilinear interplotation func
function result = lerp(t, a, b)
    result = a + t .* (b - a);
end
