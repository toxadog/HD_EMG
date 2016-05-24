function [Iu,Id] = imagecut(I,a,b)
%D?couper l'image I
%   a=[x1,y1]; b=[x2,y2];
Il=[];
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
        if j*K(1)+K(2)>=i
            Iu(i,j)=I(i,j);
        else
            Id(i,j)=I(i,j);
        end
    end
end
     
        

end

