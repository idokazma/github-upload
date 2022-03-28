        hold on 
        for k=1:length(params.sca_x)
            
        theta = 0 : (2 * pi / 10000) : (2 * pi);
        pline_x = params.radius * cos(theta) + params.sca_x(k);
        pline_y = params.radius * sin(theta) + params.sca_y(k);
        plot(pline_x, pline_y, '-r');
        
        theta = 0 : (2 * pi / 10000) : (2 * pi);
        pline_x = params.radius*1.2 * cos(theta) + params.sca_x(k);
        pline_y = params.radius*1.2 * sin(theta) + params.sca_y(k);
        plot(pline_x, pline_y, '-r');
        
        theta = 0 : (2 * pi / 10000) : (2 * pi);
        pline_x = params.radius*0.8 * cos(theta) + params.sca_x(k);
        pline_y = params.radius*0.8 * sin(theta) + params.sca_y(k);
        plot(pline_x, pline_y, '-r');
        end
        hold off