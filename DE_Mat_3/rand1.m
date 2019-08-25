function [V]=rand1(X,M,F,m)

	R=randperm(M);
	j=R(1);
	k=R(2);
	p=R(3);
	u=R(4);
	v=R(5);
	if j==m
	   j=R(6);
	elseif k==m
	   k=R(6);
	elseif p==m
	   p=R(6);	
	elseif u==m
	   u=R(6);	
	elseif v==m
	   v=R(6);					   
	end
	V=X(j,:)+F*(X(k,:)-X(p,:));