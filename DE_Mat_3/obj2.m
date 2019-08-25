function f=obj2(X)
Nk=1;
hik=1;
dik=.251;
k=6;
%m=0;
r=50;
c=2;
%y=[150 90; -100 -120; -80 130; 140 -70;
    %60 120; -90 -130];
y=[85,95];

for j=1:r
    for i=1:c
       m(j,i)=m+((1/hik)^2*(sqrt((95-X(j))^2+(95-X(i)^2))-dik)^2);
    end
end
f=(1/Nk)*m;
end


