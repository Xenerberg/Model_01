q_ref = [0;0;0;1];
q_real = [0.0200;0.0200;0.0200;0.99];
q_ref_star = [-q_ref(1:3);q_ref(4)];

del_q = fn_CrossTensor(q_real,0)*q_ref_star;
