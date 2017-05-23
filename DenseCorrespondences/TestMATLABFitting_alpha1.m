function [f,p] = TestMATLABFitting_alpha1(x,y)
%{
DATE: 2012 12 30
ACA_Fitting_Sub_bata1
ACA_Fitting_bata1的子程序
最小二乘法求曲线

DATE: 2013 1 2
ACA_Fitting_Sub_bata2
ACA_Fitting_bata2的子程序
3次拟合
%}
%{
x = [0;5;10;15;20;23.15];
y = [0;-0.5;-1.5;-2.41;-0.2;0];

%}

%x=[1;1.5;2;2.5;3];y=[0.9;1.7;2.2;2.6;3];
%x = [0;5;15.25;23.15];y=[0;-0.5;-2.41;0];
%p=fittype('a*x.^2+b*x.^4+c*x.^6','independent','x');
p=fittype('a*x.^2+b*x.^4','independent','x');
f=fit(x,y,p);
%figure,plot(f,x,y); 