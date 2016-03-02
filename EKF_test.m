function [] = fn_main()
    close all;
    options = odeset('RelTol',1e-4,'AbsTol',1e-4);
    x_pre = zeros(101,1);
    x_post = zeros(101,1);
    sigma_x_pre = 0.1;
    delta_t = 0.1;
    init = 1;
    x_pre(1) = init;
    x_post(1) = init;
    Q = 0.01;
    for i_Count_1 = 2:100        
        [T,Y] = ode45(@fn_diff,[(i_Count_1-1)*delta_t,i_Count_1*delta_t],init,options);
        x_pre(i_Count_1) = init + Y(end);
        phi = cos(init);
        sigma_x_pre(i_Count_1) = phi*sigma_x_post(i_Count_1)*phi' + Q;
        init = x_pre(i_Count_1);
    end
    
    %[T,Y] = ode45(@fn_diff,[0,300],0.1,options);
    %plot(Y);
end

function dy = fn_diff(t,y)
    dy = sin(y);
end