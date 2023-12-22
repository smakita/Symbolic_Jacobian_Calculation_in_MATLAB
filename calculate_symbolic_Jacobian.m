clear
% select standard or modified DH parameters
standardDH = false;
modifiedDH = true;
n_joint = 7;
dim = 3;
% DH parameters of the manipulator when we calculate the numeric Jacobian matrix 
% dha=[0;-0.24365;-0.21325;0;0;0];
% dhalpha=[pi/2;0;0;pi/2;-pi/2;0];
% dhd=[0.1519;0;0;0.11235;0.08535;0.0819];
% dhth=[-0.209563211,-1.314294034,1.644903294,-0.341347022,1.361587612,2.889331506];

% 同次変換行列の初期化 / Initialization of homogeneous transformation matrix
alpha = sym('ap',[1 n_joint]); %alpha
theta = sym('th',[1 n_joint]); %theta
a = sym('A',[1 n_joint]);
d = sym('D',[1 n_joint]); 
zero = sym('0');
one = sym('1');

Ttotal = [one,zero,zero,zero;zero,one,zero,zero;zero,zero,one,zero;zero,zero,zero,one];
for i=1:n_joint
    % パラメータを代入する場合．文字列で計算するときはコメントにする．
%    a(i) = dha(i);
%    alpha(i) = dhalpha(i);
%    d(i) = dhd(i);
%    theta(i) = dth(i);
    Ta = [one,zero,zero,a(i);zero,one,zero,zero;zero,zero,one,zero;zero,zero,zero,one];
    Talpha=[one,zero,zero,zero;zero,cos(alpha(i)),-sin(alpha(i)),zero;zero,sin(alpha(i)),cos(alpha(i)),zero;zero,zero,zero,one];
    Td = [one,zero,zero,zero;zero,one,zero,zero;zero,zero,one,d(i);zero,zero,zero,one];
    Ttheta=[cos(theta(i)),-sin(theta(i)),zero,zero;sin(theta(i)),cos(theta(i)),zero,zero;zero,zero,one,zero;zero,zero,zero,one];
    if(standardDH == true)
        Ttotal = Ttotal*Td*Ttheta*Ta*Talpha;
    else(modifiedDH == true)
        Ttotal = Ttotal*Talpha*Ta*Ttheta*Td;
    end

end
L = [zero;zero;zero;one] % position of the end tool of the manipulator
Ttotal = simplify(Ttotal)
L = Ttotal*L;
L = [L(1);L(2);L(3)]
%Lval = subs(L,theta,dhth) %ベクトルL内のthetaに dhthの数値を代入
%pos = [double(Lval(1));double(Lval(2));double(Lval(3))] %手先位置が正しいかどうかの確認のため

J=jacobian(L, theta) %ヤコビ行列のシンボリック計算 / Calculation of the symbolic Jacobian matrix
diary jacobian.txt % save the output to a text file
for i =1:dim
    for j = 1:n_joint 
        fprintf('\n jacobian[%d][%d]\n =', i, j)
        disp(J(i,j)) 
        fprintf(';')
    end
end
diary off;