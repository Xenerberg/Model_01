function [] = fn_main()
    close all;
    fprintf('%f',5);
    options = odeset('RelTol',1e-4,'AbsTol',1e-4);
    x_pre = zeros(1001,1);
    x_post = zeros(1001,1);
    sigma_x_pre = zeros(1001,1);
    sigma_x_post = zeros(1001,1);
    K = zeros(1000,1);
    delta_t = 0.01;
    init = 1;
    x_pre(1) = 3;
    x_post(1) = 3;
    sigma_x_pre(1) = 0.01;
    sigma_x_post(1) = 0.01;
    Q = 0.1; R = 0.1;
    
    t = 0.0:delta_t:10;
    %z = cos(t);
    
    noise_process = sqrt(Q)*randn(1000,1);
    
    [T,Y1] = ode45(@fn_diff,t,init,options);
    %plot(Y1);
    z_noised = cos(Y1) + sqrt(R)*randn(1001,1);
    
    for i_Count_1 = 2:1001        
         %Prediction stage
        [T,Y] = ode45(@fn_diff,[(i_Count_1-2)*delta_t,(i_Count_1-1)*delta_t],x_post(i_Count_1 - 1),options);
        x_pre(i_Count_1) = x_post(i_Count_1 - 1) + ((Y(end)-x_post(i_Count_1-1))); %+noise_process(i_Count_1-1));
        %x_pre(i_Count_1) = Y(end);
        phi = cos(x_post(i_Count_1 - 1));
        sigma_x_pre(i_Count_1) = phi*sigma_x_post(i_Count_1-1)*phi' + Q;
        
        %Measurement stage
        residual = z_noised(i_Count_1) - cos(x_pre(i_Count_1));
        H = -sin(x_pre(i_Count_1));
        S = H*sigma_x_pre(i_Count_1)*H' + R;
        K(i_Count_1-1) = sigma_x_pre(i_Count_1)*H'*inv(S);
        x_post(i_Count_1) = x_pre(i_Count_1) + K(i_Count_1-1)*residual;
        sigma_x_post(i_Count_1) = (1 - K(i_Count_1-1)*H)*sigma_x_pre(i_Count_1);
    end
    x_pre = mod(x_pre,2*pi);
    x_post = mod(x_post,2*pi);
    %[T,Y] = ode45(@fn_diff,[0,300],0.1,options);
    %plot(Y);
    %Plotting tools
    figure;
    plot(Y1,'linewidth',3);hold all;
    plot(x_post);
    grid on;
    legend('True state','Final estimated state');
    figure;
    grid on;
    plot(x_pre);hold all;
    plot(x_post);
    legend('x_pre','x_post');
    figure;
    grid on;
    plot(sigma_x_pre);hold all;
    plot(sigma_x_post);
    legend('sigma_x_pre','sigma_x_post');
    figure;
    plot(K);
    legend('Kalman gain');
    
end

function dy = fn_diff(t,y)
    
    dy = sin(y);
end