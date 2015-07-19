function DigitCubeSum()
%Obtains all 3-digit integers that are equal to the sum of the cubes of
%their digits.

for i=1:9
    for j=0:9
        for k=0:9
            n=100*i+10*j+k;
            if(n==i^3+j^3+k^3)
                disp([num2str(n),' = ',num2str(i),'^3+',num2str(j),'^3+',num2str(k),'^3']);
            end
        end
    end
end

