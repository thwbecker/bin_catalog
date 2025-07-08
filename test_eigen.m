a=ones(3,3);
a(1,1)=1;
a(1,2)=a(2,1)=2;
a(1,3)=a(3,1)=3;
a(2,2)=4;
a(2,3)=a(3,2)=5;
a(3,3)=6;

[v d ] = eigs(a);

ev = diag(d);
[evs ind]=sort(ev);

for i=1:3
  j=ind(i);
  out = [ ev(j) v(j,:) ]
end 
