function is_sucess = plot_shadings(shadings, time, shadings_height, color_1, color_2)

face_alpha = 0.1;

if length(shadings) == length(time)
    % Mark intervals of contraction
    if color_1 ~= 0
        j_change = 1;
        for j = 1:length(time)
            if shadings(j_change) ~= shadings(j)
                x_box = [time(j_change), time(j_change), time(j), time(j)];
                y_box = [0, shadings_height, shadings_height, 0];
                if shadings(j) == 0
                    patch(x_box', y_box', 1,'FaceAlpha', face_alpha, 'EdgeColor', 'none', 'Facecolor', color_1)
                end
                j_change = j;
            end
        end
        
        x_box = [time(j_change), time(j_change), time(end), time(end)];
        y_box = [0, shadings_height, shadings_height, 0];

        if shadings(end) == 1
            patch(x_box', y_box', 2,'FaceAlpha', face_alpha, 'EdgeColor', 'none', 'Facecolor', color_1)
        end
    end
    
    if color_2 ~= 0
        j_change = 1;
        for j = 1:length(time)
            if shadings(j_change) ~= shadings(j)
                x_box = [time(j_change), time(j_change), time(j), time(j)];
                y_box = [0, shadings_height, shadings_height, 0];
                if shadings(j) == 1
                    patch(x_box', y_box', 1,'FaceAlpha', face_alpha, 'EdgeColor', 'none', 'Facecolor', color_2)
                end
                j_change = j;
            end
        end

        x_box = [time(j_change), time(j_change), time(end), time(end)];
        y_box = [0, shadings_height, shadings_height, 0];

        if shadings(end) == 0
            patch(x_box', y_box', 2,'FaceAlpha', face_alpha, 'EdgeColor', 'none', 'Facecolor', color_2)
        end
    end
    is_sucess = true;
else
    is_sucess = false;
end

end

