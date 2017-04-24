function [Nfinal,upperBoundOut,upperBoundOut2] = computeNForFixedBeta(beta,Nsystems,epsilon)
close all;
NupperBound = 7000;



tic;
%epsilon > 1
for i = 1:Nsystems
    counter2 = 1;
    for n = 2:2:6   
        counter1 = 1;
       for m = 2:2:6
        A = cell(1,m); %cell A stores the modes of the system
        d = n*(n+1)/2+1;
        
        %number of points in the sampling, between d+1 and 1000
        if (m==2 && n==2)
           N = ceil(d*1.7^(m+n));
        else
           N = 30;
        end
        X = cell(1,N);
        Y = cell(1,N);
        
        %generate the m modes
        for j = 1:m
            v=-2+rand(1,n)*2;
            D = diag(v); % random eigenvalues of the mode i
            V = orth(randn(n));
            E = V*D*V';
            A{j} = E;
        end
        
        for j= 1:NupperBound 
            v=randn(n,1);
            X{j}=v/sqrt(sum(v.^2));
            k=unidrnd(m); %random uniform generation of the index of the mode applied to sampled point X{j}
            Y{j}=A{k}*X{j};
        end
        
        
       
        
        c=jsr_prod_bruteForce(A); %computation of the JSR with brute force method provided by the JSR Toolbox
        upperBoundOut{i}(counter1, counter2) = c(2);
        upperBound = inf;
        flag = 0;
        %this should be sequential
         while (upperBound/c(2)>epsilon && flag == 0)
               [lowerBound, upperBound] =  computeRhoBlackboxGivenXY(A, beta, X(1:min(N,NupperBound)), Y(1:min(N, NupperBound)))
               upperBoundOut2{i}(counter1,counter2) = upperBound;
               N
               m
               n
               N = ceil(N*1.25);
              if (N>7000*1.25)
                  flag = 1;
              end
         end
          Nfinal{i}(counter1,counter2)=ceil(N/1.25)
          if (upperBound/c(2)>epsilon)
             Nfinal{i}(counter1,counter2)=inf;
          end
          counter1 = counter1+1;
       end
       counter2 = counter2+1;
    end    
end
toc;




end