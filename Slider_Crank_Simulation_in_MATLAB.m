clear all; 
close all;
clc;

while true
    fprintf('...Welcome to the SLIDERCRANK Matlab Simulation\nby the Team "4 Guys Solving Assignment Questions"...\n\n')
    fprintf('Please choose an option below:\n');
    fprintf('1. Do you want to enter your own input?\n');
    fprintf('2. Do you want to test the default inputs?\n');
    fprintf('3. Exit\n');

    choice1 = input('Please enter the option number: ');

    if choice1 == 1
        a = input('Enter the length of crank (in m) : ');
        b = input(['-----------------\n' ...
            'Enter the length of connecting rod (in m) : ']);
        h = input(['----------------\n' ...
            'Enter slider offset distance from crank origin (if not then enter 0) : ']);
        theta = input(['-----------------\n' ...
            'Enter the MAXIMUM ANGLE (in degrees) until which the crank should move (in anticlockwise direction)\n' ...
            '[Please enter in multiples of 360 if you want complete rotation of the Crank (For example: 360, 720, 1080, 1440, etc...)]: ']);
        k = input(['----------------\n' ...
            'Enter angular velocity of the crank (in rad/s) : ']);

    elseif choice1 == 2
        fprintf(['--------------------------------\nHere a, b, h, theta, and k are described as follows.\na = Length of the Crank' ...
            '\nb = Length of the Connecting Rod\nh = Distance of the Offset\ntheta = The Angle (in degrees) by which the crank should move ' ...
            '(in anticlockwise direction)\nk = Angular Velocity of the Crank\n----------------------------------\n'])
        fprintf('So, please choose any one of the following default values:\n');
        fprintf('1. a = 10 , b = 40, h = 0, theta = 360, k = 1\n    (This case shows the NORMAL/REGULAR working of the simulation)\n');
        fprintf('2. a = 40, b = 10, h = 0, theta = 360, k = 1\n    (This case shows the simulation when the Crank length > Connecting Rod)\n');
        fprintf('3. a = 10, b = 40, h = 5, theta = 360, k = 1\n    (This case shows the simulation when OFFSET exists)\n');
        fprintf('4. a = 10, b = 10, h = 0, theta = 360, k = 1\n    (This case shows the simulation when a = b)\n');
        fprintf('5. a = 10, b = 40, h = 0, theta = 120, k = 1\n    (This case shows the simulation for 120 degrees angle)\n');
        fprintf("6. a = 1000, b = 1420, h = 0, theta = 360, k = 1\n    (This case shows the simulation's capabilty of adjusting the axes)\n");
        fprintf("7. a = 10, b = 40, h = 0, theta = 360, k = 3\n    (This case shows how the velocity and acceleration values of the slider " + ...
            "varies with angular velocity)\n");
        fprintf('8. Exit\n');
    
        choice2 = input('Enter the option number: ');

        if choice2 == 1
            a = 10; b = 40; h = 0; theta = 360; k = 1;

        elseif choice2 == 2
            a = 40; b = 10; h = 0; theta = 360; k = 1;

        elseif choice2 == 3
            a = 10; b = 40; h = 5; theta = 360; k = 1;

        elseif choice2 == 4
            a = 10; b = 10; h = 0; theta = 360; k = 1;

        elseif choice2 == 5
            a = 10; b = 40; h = 1; theta = 120; k = 1;

        elseif choice2 == 6
            a = 1000; b = 1420; h = 0; theta = 360; k = 1;

        elseif choice2 == 7
            a = 10; b = 40; h = 0; theta = 360; k = 3;

        elseif choice2 == 8
            fprintf('---------------\nTHANK YOU\n-From Team "4 Guys Solving Assignment Questions" - Prudhvi Nallagatla, G.O.V Umesh Chandra, Pranav Sudheer, R. Sai Datta Praneeth\n--------------------\n')
            break;
        else
            fprintf('Invalid option. Please choose a valid option.\n');
            continue;
        end

    elseif choice1 == 3
        fprintf('---------------\nTHANK YOU\n-From Team "4 Guys Solving Assignment Questions" - Prudhvi Nallagatla, G.O.V Umesh Chandra, Pranav Sudheer, R. Sai Datta Praneeth\n--------------------\n')
        break;

    else
        fprintf('Invalid option. Please choose a valid option.\n');
        continue;
    end
    
    p1 = [0 0]; 
    n = b/a;
    T = linspace(0, theta, theta + 1);

    stopSimulation = false;

    if n <= 1
        fprintf("--------------------------------\nThe connecting rod's length is either SMALLER or EQUAL to the Crank's length." + ...
            "\nSo, this mechanism will not work. \nPlease try different values.\n----------------------------------------\n")
        stopSimulation = true;
    end
    mass_crank = 1;
    mass_connecting_rod = 1;
    mass_slider = 1;

    if ~stopSimulation
        for i = 1:length(T)
            
            p2 = a * [cosd(T(i)) sind(T(i))];
            p3 = [p1(1) + p2(1) + sqrt(b^2 - (p2(2))^2), h];
            
            term1 = 1 - cosd(T(i));
            term2 = (n - sqrt(n^2 - sind(T(i))^2));
            
            position(i) = a * (term1 + term2);
            velocity(i) = a * k * (sind(T(i)) + sind(2 * T(i)) / (2 * sqrt(n^2 - sind(T(i))^2)));
            %fp = w2 r(cos θ + cos2θ/n) acc of slider formula
            %wpc = w cos θ / (n2 – sin2θ)1/2 ang vel of connecting rod
            %αpc =  -w2sinθ (n2 -1)/ (n2 -sin2θ)3/2 ang acc of connecting rod
            accl(i) = a * k^2 * (cosd(T(i)) + cosd(2*T(i))/n);
            angvel_cr(i) = (a*k*cosd(T(i)))/sqrt(n*n - sind(2*T(i)));
            angacc_cr(i) = (-a*k*k*sind(T(i)*(n*n - 1)))/(n*n-sind(2*T(i)))^3/2;

            
            figure(1)
            set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1])
            plot([p1(1) p2(1)], [p1(2) p2(2)], 'linewidth', 10);
            axis([-a - 2, a + b + 2, -(a + b + 2), a + b / 2 + 2]);
            hold on
            
            ang = linspace(0, 2 * pi);
            plot(a * cos(ang), a * sin(ang), 'k--')
           
            plot([p2(1) p3(1)], [p2(2) p3(2)], 'linewidth', 10); 
            plot(p1(1), p1(2), 'marker', 'o', 'markersize', 6, 'markerfacecolor', 'k');
            plot(p2(1), p2(2), 'marker', 'o', 'markersize', 6, 'markerfacecolor', 'k');
            plot(p3(1), p3(2), 'marker', 's', 'markersize', 30, 'markerfacecolor', 'g');
            % Velocity analysis
            % Position points
            p1 = [0 0];
            crank_end = a * [cosd(T(i)) sind(T(i))];
            slider_end = [position(i) h];


            % Velocity vectors with arrowheads
            arrowhead_size = 0.75;  % Adjust the arrowhead size as needed

            % Calculate unit vectors for the direction of motion
            crank_dir = (crank_end - p1) / norm(crank_end - p1);
            slider_dir = (slider_end - crank_end) / norm(slider_end - crank_end);
            perpendicular_dir = [crank_dir(2)-slider_dir(2), crank_dir(1) +slider_dir(1)];
            connecting_rod_dir = perpendicular_dir;
            

            % Calculate the perpendicular velocity vectors
            crank_vel = a * k * [-crank_dir(2), crank_dir(1)];
            slider_vel = velocity(i) * [-1, 0];
            connecting_rod_vel= slider_vel - crank_vel;
            slider_start = (a+b) - position(i);  % Calculate the starting point for the red vector

            % Create velocity vectors with lengths and arrowheads
            %quiver(slider_start, h , slider_vel(1)/k, slider_vel(2)/k, 'r', 'LineWidth', 2, 'MaxHeadSize', arrowhead_size, 'AutoScalefactor', 1);
            %quiver(crank_end(1), crank_end(2), connecting_rod_vel(1)/k, connecting_rod_vel(2)/k, 'color', "#EDB120", 'LineWidth', 2, 'MaxHeadSize', arrowhead_size, 'AutoScalefactor', 1);
            %quiver(crank_end(1), crank_end(2), crank_vel(1)/k, crank_vel(2)/k, 'b', 'LineWidth', 2, 'MaxHeadSize', arrowhead_size, 'AutoScalefactor', 1);

            %velocity diagram
            text(a - 0.1*b, a - 0.65*b, "Velocity Vector Polygon:", 'EdgeColor', 'r') 
            origin_point = [ a + 2, -(a + b/2) ];
            quiver(origin_point(1), origin_point(2), slider_vel(1)/k, slider_vel(2)/k, 'r', 'LineWidth', 2, 'MaxHeadSize', arrowhead_size, 'AutoScalefactor', 1);
            quiver(origin_point(1), origin_point(2), crank_vel(1)/k, crank_vel(2)/k, 'b', 'LineWidth', 2, 'MaxHeadSize', arrowhead_size, 'AutoScalefactor', 1);

            quiver(origin_point(1)+ crank_vel(1)/k, origin_point(2) + crank_vel(2)/k, connecting_rod_vel(1)/k, connecting_rod_vel(2)/k,  'color', "#EDB120", 'LineWidth', 2, 'MaxHeadSize', arrowhead_size, 'AutoScalefactor', 1);
        
            
            % Acceleration analysis
            vel_crank_mag = k * a;
            acc_crank_mag = (vel_crank_mag^2) / a;
            acc_crank = [-acc_crank_mag * cosd(T(i)), -acc_crank_mag * sind(T(i))];
            
            connecting_rod_vel_mag = sqrt(connecting_rod_vel(1)^2 + connecting_rod_vel(2)^2);
            acc_crod_centri_mag = (connecting_rod_vel_mag^2) / b;
            %quiver(p2(1),p2(2),acc_crank(1)/(k*k),acc_crank(2)/(k*k), "Color", 'k','LineWidth', 2);
            

            % acc plot
            text(a + 0.6*b, a - 0.65*b, "Acceleration Vector Polygon:", 'EdgeColor', 'r') 
            op = [ a + 1.1*b, a - b];
            quiver(op(1),op(2),acc_crank(1)/(k*k),acc_crank(2)/(k*k),'k','LineWidth', 2,'AutoScalefactor', 1);
            
            
            beta = asind(sind(T(i))/n);
            acc_crod_centri = [-1* acc_crod_centri_mag * cosd(beta), acc_crod_centri_mag * sind(beta)];
            acc_crod_start = [op(1)+(acc_crank(1)/(k*k)), op(2)+(acc_crank(2)/(k*k))];
            quiver(acc_crod_start(1), acc_crod_start(2), acc_crod_centri(1)/(k*k), acc_crod_centri(2)/(k*k),"Color", "#EDB120", 'LineWidth', 2, 'AutoScalefactor', 1, 'LineStyle', ':', "MaxHeadSize", arrowhead_size);
            
            acc_slider = accl(i) * [-1, 0];
            quiver(op(1), op(2), acc_slider(1)/(k*k), acc_slider(2)/(k*k),"Color", 'r','LineWidth', 2, 'AutoScalefactor', 1, "MaxHeadSize", arrowhead_size);
            %quiver(slider_start, h, acc_slider(1)/(k*k), acc_slider(2)/(k*k),"Color", 'r','LineWidth', 2, 'AutoScalefactor', 1, "MaxHeadSize", arrowhead_size);
            
            acc_crod_tang = acc_slider - acc_crank - acc_crod_centri;
            acc_crod_tang_start = acc_crod_start + (acc_crod_centri/(k*k));
            acc_crod_tang_mag = sqrt(acc_crod_tang(1)^2 + acc_crod_tang(2)^2);
            quiver(acc_crod_tang_start(1), acc_crod_tang_start(2), acc_crod_tang(1)/(k*k), acc_crod_tang(2)/(k*k), "Color", "#0072BD", 'LineWidth', 2, 'AutoScalefactor', 1, 'LineStyle', ':', "MaxHeadSize", arrowhead_size);
            
            acc_crod_mag = sqrt(acc_crod_centri_mag^2 + acc_crod_tang_mag^2);
            acc_crod = acc_crod_tang+acc_crod_centri;
            quiver(acc_crod_start(1), acc_crod_start(2), acc_crod(1)/(k*k), acc_crod(2)/(k*k), "Color", "#77AC30", 'LineWidth', 2,'AutoScalefactor', 1, "MaxHeadSize", arrowhead_size);
            %quiver(p2(1),p2(2), acc_crod(1)/(k*k), acc_crod(2)/(k*k), "Color", "#77AC30", 'LineWidth', 2,'AutoScalefactor', 1, "MaxHeadSize", arrowhead_size);

            %Static force analysis:
            % Initialize force vectors
            force_crank_tang = zeros(size(T));
            force_crank_norm = zeros(size(T));
            force_crod_tang = zeros(size(T));
            force_crod_norm = zeros(size(T));
            force_slider = zeros(length(T), 2);

            % Calculate acceleration components for force vector calculation
            acc_crank_tang = -a * k^2 * sind(T);  % Tangential acceleration of the crank
            acc_crank_norm = a * k^2 * (1 - cosd(T));  % Normal acceleration of the crank
            
            acc_crod_tang = -a * k^2 * sind(T);  % Tangential acceleration of the connecting rod
            acc_crod_norm = a * k^2 * (cosd(T));  % Normal acceleration of the connecting rod
            
            % Calculate force vectors in the resultant direction of normal and tangential components
            force_crank_tang = -mass_crank * acc_crank_tang;  % Negative sign to align with acceleration
            force_crank_norm = -mass_crank * acc_crank_norm;  % Negative sign to align with acceleration
            force_crank_resultant = force_crank_tang + force_crank_norm;
            
            force_crod_tang = -mass_connecting_rod * acc_crod_tang;  % Negative sign to align with acceleration
            force_crod_norm = -mass_connecting_rod * acc_crod_norm;  % Negative sign to align with acceleration
            force_crod_resultant(:, 1) = -mass_connecting_rod * acc_crod_tang .* cosd(T);
            force_crod_resultant(:, 2) = -mass_connecting_rod * acc_crod_norm .* sind(T);
            
            % Calculate force vector for the slider in the direction of its acceleration
            force_slider = -mass_slider * accl;
            
            % Plot force vectors
            midpoint_crank = 0.5 * (p1 + p2);
            midpoint_crod = 0.5 * (p2 + p3);
            midpoint_slider = [slider_start, h];

            % Crank force vector
            quiver(midpoint_crank(1), midpoint_crank(2), force_crank_resultant(i), 0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', arrowhead_size, 'AutoScalefactor', 1);

            % Connecting Rod force vector
            quiver(midpoint_crod(1), midpoint_crod(2), force_crod_resultant(i, 1), force_crod_resultant(i, 2), 'Color', 'g', 'LineWidth', 2, 'MaxHeadSize', arrowhead_size, 'AutoScalefactor', 1);
            
            % Slider force vector
            quiver(midpoint_slider(1), midpoint_slider(2), force_slider(i), 0, 'Color', 'b', 'LineWidth', 2, 'MaxHeadSize', arrowhead_size, 'AutoScalefactor', 1);

            %Plot table
            val2 = [sprintf('Position Analysis: '), ...
               sprintf(['-------------------------------------------------\n' ...
               'Crank Length: %5.2fm'], a), ... 
               sprintf('\nConnecting Rod Length: %5.2fm', b), ...
               sprintf('\nOffset Distance: %5.2fm', h), ...
               sprintf('\nCrank angle (in degrees): %d', round(T(i))), ...
               sprintf("\nSlider's Instantaneous Position: %5.3f", (a+b)-position(:,i)), ...
               sprintf("\nSlider's Instantaneous Velocity (in m/s): %5.3f", (velocity(:,i)))];
            
            fontsize(10, "points")
            text(a - 1.4*b, a, val2, 'EdgeColor', 'b') 

            %Plot table 2 - This table is for the velcity vector diagram
            %details
            val3 = [sprintf('Velocity Vector Analysis: '), ...
               sprintf("-----------------------------------------------\n" + ...
               "Connecting Rod's Angular Velocity: %5.2f rad/sec", angvel_cr(:,i)), ... 
               sprintf("\nSlider's Instantaneous Velocity: %5.3f m/sec", (velocity(:,i))), ...
               sprintf("\nGiven Crank's Angular Velocity: %5.2f rad/sec", k), ...
               sprintf('\n-----------------------------------------------\nColour of the Crank is BLUE'), ...
               sprintf('\nColour of the Connecting Rod is YELLOW'), ...
               sprintf('\nColour of the Slider is GREEN')];
            
            fontsize(10, "points")
            text(a - 1.4*b, a - b, val3, 'EdgeColor', 'b')


            %Plot table 3 - This table is for the Acceleration vector diagram
            val1 = [sprintf('Acceleration Vector Analysis: '), ...
               sprintf("-----------------------------------------------\nSlider's Instantaneous Acceleration: %5.3f m/s^2", accl(:,i)), ... 
               sprintf("\nCrank's Angular Acceleration: %5.1f rad/s^2", 0), ...
               sprintf("\nConnecting Rod's Angular Acceleration: %5.5f rad/s^2", angacc_cr(:, i)), ...
               sprintf("\n-----------------------------------------------\n" + ...
               "Colour of the Crank's\n      Normal Accleration component is BLACK"), ...
               sprintf("\nColour of the Connecting Rod's\n      Centripetal Accleration component is YELLOW"), ...
               sprintf("\nColour of the Connecting Rod's\n      Tangential Accleration component is BLUE"), ...
               sprintf("\nColour of the Connecting Rod's\n      Resultant Accleration is GREEN"), ...
               sprintf("\nColour of the Slider's Accleration is RED")];
            
            fontsize(10, "points")
            text(a + 1.4*b, a-0.25*b, val1, 'EdgeColor', 'b')
            
            
            axis equal
            hold off
        end

        %---------------Table's Code----------------------
        poscol = zeros(size(T));
        velcol = zeros(size(T));
        acclcol = zeros(size(T));
        cr_angvelcol = zeros(size(T));
        cr_angacclcol = zeros(size(T));
        
        for i = 1:numel(T)
            poscol(i) = position(i);
            velcol(i) = velocity(i);
            acclcol(i) = accl(i);
            cr_angvelcol(i) = angvel_cr(i);
            cr_angacclcol(i) = angacc_cr(i);
        end
        
        data = [round(T'), poscol', velcol', acclcol', cr_angvelcol', cr_angacclcol'];
        data = num2cell(data);
        h = {'Crank Angle', "Slider's Position", "Slider's Velocity", "Slider's Accleration", ...
            "Coupler's Angular Velocity", "Coupler's Angular Accleration"}.'; 
        
        f = figure;  
        columnWidths = {'auto', 'auto', 'auto', 'auto', 'auto', 'auto'};
        t = uitable(f, 'Data', data, 'ColumnName', h, 'Units','normalized','Position', [0 0 1 1]);

        
    end

    choice3 = input(['Thank you for simulating SliderCrank Mechanism.\nNow, what do you want to do?\n' ...
        '1. Do you want to try the simulation again? (OR)\n' ...
        '2. Do you want to EXIT the program?\n' ...
        'Please enter the option number (NOTE: "2" will STOP the simulation): ']);

    if choice3 == 2
        fprintf('---------------\nTHANK YOU\n-From Team "4 Guys Solving Assignment Questions" - Prudhvi Nallagatla, G.O.V Umesh Chandra, Pranav Sudheer, R. Sai Datta Praneeth\n--------------------\n')
        break;
    end
end
