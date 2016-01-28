%Function to compute Matrix exponential
function [exp_Matrix] = fn_CreateMatrixExp(M)
    n = length(M);
    switch(n)
        case 3       
            v_eig = eig(M);
            M_eig = [ones(3,1),v_eig,v_eig.^2];
            v_rhs = exp(v_eig);

            v_m = (M_eig)\v_rhs;
            exp_Matrix = v_m(1)*eye(3,3) + v_m(2)*M + v_m(3)*M*M;
        case 4
            v_eig = eig(M);
            M_eig = [ones(4,1),v_eig,v_eig.^2,v_eig.^3];
            v_rhs = exp(v_eig);

            v_m = (M_eig)\v_rhs;
            exp_Matrix = v_m(1)*eye(4,4) + v_m(2)*M + v_m(3)*M*M + v_m(4)*M*M*M;
        otherwise
            exp_Matrix = 0;
    end
end