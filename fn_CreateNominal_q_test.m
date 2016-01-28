function q_nominal = fn_CreateNominal_q_test(omega,t_delta,n,q)
%#codegen
    omega_norm = norm(omega);
    c = cos(omega_norm*t_delta/2);
    s = sin(omega_norm*t_delta/2);
    exponential_first =( c + s )*eye(4,4);
    q_omega = [omega;0];
    n_v = [0;0;n];
    q_n = [n_v;0];
    exponential_second = ((omega_norm*t_delta*c - 4*s)/(2*omega_norm))*fn_CrossTensor(q_omega,0);
    exponential_om = exponential_first - exponential_second;
    exponential_n = eye(4,4);%exp((1/2)*t_delta*(-fn_CrossTensor(q_n,1)));
    exponential = exponential_om*exponential_n;
    q_nominal = exponential*q;
end