function [deltaBoundFrac, deltaBound] = computeGapForFixedBetaAndN(beta,Nsystems,N)
close all;

tic;
deltaBound = cell(1,Nsystems);
deltaBoundFrac = cell(1, Nsystems);
%epsilon > 1
for i = 1:Nsystems
    counter2 = 1;
    for n = 2:2:6   
        counter1 = 1;
       for m = 2:2:6
        A = cell(1,m); %cell A stores the modes of the system
        d = n*(n+1)/2+1;
        
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
        
        for j= 1:N
            v=randn(n,1);
            X{j}=v/sqrt(sum(v.^2));
            k=unidrnd(m); %random uniform generation of the index of the mode applied to sampled point X{j}
            Y{j}=A{k}*X{j};
        end
      
        c=jsr_prod_bruteForce(A); %computation of the JSR with brute force method provided by the JSR Toolbox
        
        %this should be sequential
        [~, upperBound] =  computeRhoBlackboxGivenXY(A, beta, X, Y);
        deltaBound{i}(counter1, counter2) = upperBound - c(2);
        deltaBoundFrac{i}(counter1, counter2) = upperBound/c(2);
        counter1 = counter1+1;
       end
       counter2 = counter2+1;
    end    
end
toc;


end