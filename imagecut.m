function [Iu,Id] = imagecut(I,a,b)
%D?couper l'image I
%   a=[x1,y1]; b=[x2,y2];
Il=[];
h=tukeywin(32,0.25);
for i=1:32
    Il=[Il I(i,:)];
end
backgr=quantile(Il,0.2);
X=[a(2) 1; b(2) 1];
Y=[a(1);b(1)];
K=X\Y;
Iu=ones(size(I)).*backgr;
Id=ones(size(I)).*backgr;
for i=1:size(I,1)
    for j=1:size(I,2)
        d=abs(-K(1)*j+i*1-K(2))/sqrt(K(1)^2+1^2);
        if j*K(1)+K(2)>=i
            Iu(i,j)=I(i,j)*h(fix(d)+1)+backgr;
        else
            Id(i,j)=I(i,j)*h(fix(d)+1)+backgr;
        end
    end
end
     
        

end

