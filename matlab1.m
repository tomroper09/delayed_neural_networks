function lowert2_eigenvalues_t1_varies()

    % Define the parameters
    alpha = 1;
    beta = 1;
    a = 10;
    b = -10;
    c = 10;
    d = 2;
    tau2 = 0.1;
    u_star = 0.1864;
    v_star = 0.1338;

    % Initialise vectors
    tau1_values = linspace(0, 12, 1000);
    num_eigenvalues = 2; % Number of eigenvalues to solve for
    num_runs = 3;

    % Initialize a matrix to store the real parts of eigenvalues
    real_parts_all_runs = zeros(length(tau1_values), num_eigenvalues * num_runs);

    options = optimoptions('fsolve', 'Display', 'none');

    % Initial guesses array
    initial_guesses = [0.4+1i, 0.4-1i; 1+1.5i, 1-1.5i ; 0+0.5i, 0-0.5i]; % Add more initial guesses as needed

    for k = 1:size(initial_guesses, 1)
        lambda_guess = initial_guesses(k, :);

        % Choose direction for continuation (normal or reverse)
        iteration_order = 1:length(tau1_values);
        if k > 1 % Reverse order for subsequent runs
            iteration_order = fliplr(iteration_order);
        end

        for i = iteration_order
            tau1 = tau1_values(i);

            % e(lambda) matrix
            E = @(lambda) [lambda + 1 - a*beta*u_star*(1 - u_star)*exp(-lambda*tau1), -b*beta*u_star*(1 - u_star)*exp(-lambda*tau2);
                           -c*beta*v_star*(1 - v_star)*exp(-lambda*tau2), (lambda/alpha) + 1 - d*beta*v_star*(1 - v_star)*exp(-lambda*tau1)];

            % Define the determinant function
            det_E = @(lambda) det(E(lambda));

            % Solve for the eigenvalues
            for j = 1:length(lambda_guess)
                lambda_eigenvalue = fsolve(det_E, lambda_guess(j), options);
                lambda_guess(j) = lambda_eigenvalue; % Update guess for next iteration

                % Store real parts of eigenvalues
                real_parts_all_runs(i, (k-1)*num_eigenvalues + j) = real(lambda_eigenvalue);
            end
        end
    end

    % Plot the real part of the eigenvalues vs tau1
    figure(1);
    clf;
    
    for j = 1:num_eigenvalues*num_runs
        plot(tau1_values, real_parts_all_runs(:, j), 'LineWidth', 3);
        hold on;
    end
    plot(tau1_values, zeros(1,length(tau1_values)), '--k', 'LineWidth', 2)
    xlabel('Tau 1', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel('Re(\lambda)', 'FontSize', 20, 'FontWeight', 'bold');
    title('Real Part of Eigenvalues vs Tau1 (Tau2 = 0.5)');
    set(gca, 'FontSize', 20, 'FontWeight', 'bold'); % Increase font size for axes
    set(gca, 'LineWidth', 2); % Increase axis line width

    ylim([-0.5 0.5])

    grid on;
    box off;
    hold off;
    set(gcf, 'Color', 'w'); % Set background color to white

    disp(real_parts_all_runs)
end
