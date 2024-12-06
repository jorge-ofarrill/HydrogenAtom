function HydrogenAtom(n, l_init, m_init)
    % Hydrogen Atom Probability Density Visualization
    close all;
    a0 = 1; % Bohr radius
    rdom = 2 * n^2 * a0; % Spatial domain scaling
    grid_size = 100; % Resolution of the grid

    % Create 3D grid
    x = linspace(-rdom, rdom, grid_size + 1);
    y = linspace(-rdom, rdom, grid_size + 1);
    z = linspace(-rdom, rdom, grid_size + 1);
    [X, Y, Z] = meshgrid(x, y, z);

    % Initialize quantum numbers
    l = l_init;
    m = m_init;

    % Precompute and visualize
    [psi_squared, density_threshold] = computeProbabilityDensity();

    % Visualize slices
    visualizeSlices();

    % Visualize isosurface
    fig = figure('Name', 'Hydrogen Atom', 'Color', 'k', 'Position', [200, 200, 800, 600]);
    axis equal
     xlim([-rdom, rdom]), ylim([-rdom, rdom]), zlim([-rdom, rdom]);
    isosurface_handle = plotIsosurface(density_threshold);
    threshold = density_threshold;
    camlight; 
    % Add slider for threshold adjustment
    slider = uicontrol('Style', 'slider', ...
        'Min', -11, 'Max', log10(max(psi_squared(:))), ...
        'Value', log10(density_threshold), ...
        'Units', 'normalized', ...
        'Position', [0.15 0.05 0.7 0.03], ...
        'Callback', @updateIsosurface);
  addlistener(slider, 'ContinuousValueChange', @updateIsosurface);
    % Add buttons for quantum number adjustment

    uicontrol('Style', 'pushbutton', ...
        'String', '+N', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.95 0.05 0.05], ...
        'Callback', @increaseN);

    uicontrol('Style', 'pushbutton', ...
        'String', '-N', ...
        'Units', 'normalized', ...
        'Position', [0.07 0.95 0.05 0.05], ...
        'Callback', @decreaseN);

    uicontrol('Style', 'pushbutton', ...
        'String', '+L', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.85 0.05 0.05], ...
        'Callback', @increaseL);

    uicontrol('Style', 'pushbutton', ...
        'String', '+M', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.75 0.05 0.05], ...
        'Callback', @increaseM);

    uicontrol('Style', 'pushbutton', ...
        'String', '-L', ...
        'Units', 'normalized', ...
        'Position', [0.07 0.85 0.05 0.05], ...
        'Callback', @decreaseL);

    uicontrol('Style', 'pushbutton', ...
        'String', '-M', ...
        'Units', 'normalized', ...
        'Position', [0.07 0.75 0.05 0.05], ...
        'Callback', @decreaseM);

    % Slider Label
    slider_label = uicontrol('Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.15 0.01 0.7 0.03], ...
        'String', ['Threshold = ', num2str(density_threshold)], ...
        'ForegroundColor', 'w', ...
        'BackgroundColor', 'k');

    % Nested functions for computation and updates
    function [psi_squared, density_threshold] = computeProbabilityDensity()
        rho = sqrt(X.^2 + Y.^2 + Z.^2);
        phi = atan2(Y, X);
        theta = acos(Z ./ max(rho, eps)); % Avoid division by zero

        % Radial wavefunction and spherical harmonics
        R_nl = radial_wavefunction(n, l, rho, a0);
        Y_lm = spherical_harmonic(l, m, theta, phi);

        % Probability density
        psi_squared = abs(R_nl .* Y_lm).^2;
        dV = (2 * rdom / grid_size)^3; % Volume element
        psi_squared = psi_squared / sum(psi_squared(:) * dV); % Normalize

        % Initial threshold
        density_threshold = max(psi_squared(:)) * 0.1;
    end



    function h = plotIsosurface(threshold)
        h = patch(isosurface(X, Y, Z, psi_squared, threshold));
        set(h, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        title({['n=', num2str(n), ', l=', num2str(l), ', m=', num2str(m)]; ...
               ['Threshold = ', num2str(threshold)]}, 'Color', 'w');
        xlabel('x (a.u.)'); ylabel('y (a.u.)'); zlabel('z (a.u.)');
        lighting gouraud;
        axis off; set(gca, 'Color', 'k');
    end

    function updateIsosurface(src, ~)
        threshold = 10^(src.Value); % Convert from log scale
        delete(isosurface_handle); % Remove old isosurface
        isosurface_handle = plotIsosurface(threshold);
        set(slider_label, 'String', ['Threshold = ', num2str(threshold)]);
    end

function increaseN(~, ~)
        n = n + 1;

        refreshVisualization();
        visualizeSlices()
end

function decreaseN(~, ~)
        n = n - 1;
        if n < 1
            n=1;
        end
        refreshVisualization();
        visualizeSlices()
    end

    function increaseL(~, ~)
        l = l + 1;
        if l > n - 1
            l = n - 1;
        end
        refreshVisualization();
        visualizeSlices()
    end

    function increaseM(~, ~)
        m = m + 1;
        if m > l
            m = l;
        end
        refreshVisualization();
        visualizeSlices()
    end

function decreaseL(~, ~)
        l = l - 1;
        if l < 1
            l = 1;
        end
        refreshVisualization();
        visualizeSlices()
    end

    function decreaseM(~, ~)
        m = m - 1;
        if abs(m) > l
            m = -l;
        end
        refreshVisualization();
        visualizeSlices()
    end

    function refreshVisualization()
        [psi_squared, density_threshold] = computeProbabilityDensity();
        delete(isosurface_handle); % Remove old isosurface
        isosurface_handle = plotIsosurface(threshold);
    end

    function visualizeSlices()
        figure(1);
        set(gcf, 'Color', 'k', 'Position', [100, 100, 1200, 400]);

        x_index = find(abs(x) < 1e-5, 1);
        y_index = find(abs(y) < 1e-5, 1);
        z_index = find(abs(z) < 1e-5, 1);

        slice_datax = nthroot(squeeze(psi_squared(x_index, :, :)), 5);
        slice_datay = nthroot(squeeze(psi_squared(:, y_index, :)), 5);
        slice_dataz = nthroot(squeeze(psi_squared(:, :, z_index)), 5);

        subplot(131), imagesc(y, z, slice_datax');
        axis equal tight; axis off; colormap(cool); set(gca, 'YDir', 'normal'); title('Slice at x=0');
        subplot(132), imagesc(x, z, slice_datay');
        axis equal tight; axis off; colormap(cool); set(gca, 'YDir', 'normal'); title('Slice at y=0');
        subplot(133), imagesc(x, y, slice_dataz');
        axis equal tight; axis off; colormap(cool); set(gca, 'YDir', 'normal'); title('Slice at z=0');
    end
end
