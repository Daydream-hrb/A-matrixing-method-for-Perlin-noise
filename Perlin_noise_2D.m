
clear;
clc;

% Config
size_of_noise = 256;
perlin_noise = zeros(size_of_noise, size_of_noise, 'single');
for i = 2: 8
    perlin_noise = perlin_noise + perlin(size_of_noise, 2^i);
end
[X, Y] = meshgrid(1: size_of_noise, 1: size_of_noise);

% Erosion (you can process the raw noise anyway you like)
lowest = min(perlin_noise, [], 'all');
[GX, GY] = gradient(perlin_noise);
G = sqrt(GX .^ 2 + GY .^ 2);
perlin_noise = (perlin_noise - lowest) .* (1 ./ (1 + 0.2 .* G)) + lowest;
perlin_noise(perlin_noise <= 0) = -(-perlin_noise(perlin_noise <= 0)) .^ 0.8;

% Visualization
amp_factor = 10;
figure('Color', 'w', 'Name', '2D Perlin Noise');
surf(X, Y, amp_factor * perlin_noise, 'EdgeColor', 'none');
colormap jet;
xlabel('\it{x}');
ylabel('\it{y}');
zlabel('\it{z}');
colorbar;
axis equal
set(gca, 'FontName', 'Times New Roman')

function perlin_noise = perlin(size_of_noise, freq)

    % Config
    xy_step = size_of_noise / freq;
    assert(mod(size_of_noise, freq) == 0)
    Amp = 10 / freq;
    
    % Generate gradients of grid nodes
    grad_mat = Amp * randn(size_of_noise / xy_step + 1, size_of_noise / xy_step + 1, 2); % 行 × 列 × 2
    
    % Coordinates of all nodes
    [X, Y] = ndgrid(1: size_of_noise, 1: size_of_noise);
    
    % Unit increment
    dx = single(mod(X - 0.5, xy_step) / xy_step); % n x n
    dy = single(mod(Y - 0.5, xy_step) / xy_step); % n x n

    % Gradient matrix of all nodes
    Gx = single(repelem(squeeze(grad_mat(:, :, 1)), xy_step, xy_step)); % n x n
    Gy = single(repelem(squeeze(grad_mat(:, :, 2)), xy_step, xy_step)); % n x n
    
    % Gradient bias list
    Gb = single([1, xy_step + 1, xy_step * freq, xy_step * (freq + 1)]);

    % Integral (most exciting part of the matrixing method which avoids any loop structure)
    dot00 = Gx(Gb(1): Gb(3), Gb(1): Gb(3)) .* dx + Gy(Gb(1): Gb(3), Gb(1): Gb(3)) .* dy; % n x n
    dot01 = Gx(Gb(1): Gb(3), Gb(2): Gb(4)) .* dx + Gy(Gb(1): Gb(3), Gb(2): Gb(4)) .* (dy - 1); % n x n
    dot10 = Gx(Gb(2): Gb(4), Gb(1): Gb(3)) .* (dx - 1) + Gy(Gb(2): Gb(4), Gb(1): Gb(3)) .* dy; % n x n
    dot11 = Gx(Gb(2): Gb(4), Gb(2): Gb(4)) .* (dx - 1) + Gy(Gb(2): Gb(4), Gb(2): Gb(4)) .* (dy - 1); % n x n

    % Smoothing
    tx = fade(dx);
    ty = fade(dy);
    
    % Bilinear interplotation
    perlin_noise = lerp(ty, lerp(tx, dot00, dot10), lerp(tx, dot01, dot11));

end

% Smoothing func
function t_ = fade(t)
    t_ = t .* t .* t .* (t .* (t .* 6 - 15) + 10);
end

% Bilinear interplotation func
function y = lerp(t, a, b)
    y = a + t .* (b - a);
end


