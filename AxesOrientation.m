function [lambda,fi]=AxesOrientation(evecs)

a_vec=evecs(:,1);
b_vec=evecs(:,2);
c_vec=evecs(:,3);

lambda_a=atan2(a_vec(2),a_vec(1))*180/pi;
lambda_b=atan2(b_vec(2),b_vec(1))*180/pi;
lambda_c=atan2(c_vec(2),c_vec(1))*180/pi;

fi_a=atan2(a_vec(3),norm([a_vec(1) a_vec(2)]))*180/pi;
fi_b=atan2(b_vec(3),norm([b_vec(1) b_vec(2)]))*180/pi;
fi_c=atan2(c_vec(3),norm([c_vec(1) c_vec(2)]))*180/pi;

lambda=[lambda_a lambda_b lambda_c];
fi=[fi_a fi_b fi_c];