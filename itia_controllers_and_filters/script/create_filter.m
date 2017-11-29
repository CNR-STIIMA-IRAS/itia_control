clear all;clc;close all
s=tf('s');

dt=1e-3;

z=tf('z',dt);
zi=1/z;

orders=[0 0 1 2];
names={'position_filter','effort_filter','velocity_filter','acceleration_filter'};

wn=70*2*pi;
xci=0.9;
s=tf('s');
low_pass_filter_c=zpk(1/(s^2/wn^2+2*xci*s/wn+1));
%low_pass_filter=c2d(ss(low_pass_filter_c),dt,'zoh');

deriv=zpk([0],[-5*wn],5*wn);
for idx=1:4
  order=orders(idx);
  name=names{idx};
  %low_pass_filter=1;
  low_pass_filter_c=deriv^order*low_pass_filter_c;
  Fd=minreal(ss(c2d(low_pass_filter_c,dt,'zoh')));

  Fd=minreal(Fd);
  
  A=Fd.a;
  B=Fd.b;
  C=Fd.c;
  D=Fd.d;
  Baw=zeros(size(A,1),size(C,1));
  ub=+1e6*ones(size(C,1));
  lb=-1e6*ones(size(C,1));
  x0=zeros(size(A,1),1);
  
  fid=1;
  fprintf(fid,'  %s:\n',name);
  fprintf(fid,'    A:\n');
  fprintf(fid,'%s',save_matrix_to_yaml(A,4));
  fprintf(fid,'    B:\n');
  fprintf(fid,'%s',save_matrix_to_yaml(B,4));
  fprintf(fid,'    Baw:\n');
  fprintf(fid,'%s',save_matrix_to_yaml(Baw,4));
  fprintf(fid,'    C:\n');
  fprintf(fid,'%s',save_matrix_to_yaml(C,4));
  fprintf(fid,'    D:\n');
  fprintf(fid,'%s',save_matrix_to_yaml(D,4));
  fprintf(fid,'    max_output:\n');
  fprintf(fid,'%s',save_matrix_to_yaml(ub,4));
  fprintf(fid,'    min_output:\n');
  fprintf(fid,'%s',save_matrix_to_yaml(lb,4));
  fprintf(fid,'    initial_state:\n');
  fprintf(fid,'%s',save_matrix_to_yaml(x0,4));
  
  bode(Fd,tf('s')^order);
end
%%
t=(0:dt:1)';
y=2*ones(size(t));
y=t.^2;
y=t;
Dy=zeros(size(t));
x=x0;
for idx=1:length(t)
  Dy(idx)=C*x+D*y(idx);
  x=A*x+B*y(idx);
end
plot(t,Dy)