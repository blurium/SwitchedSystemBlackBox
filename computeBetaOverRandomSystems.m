function [deltaC,proba] = computeBetaOverRandomSystems(beta,Nsystems,mUpperBound,nUpperBound,NupperBound)
counterOK = 0; 
for i = 1:Nsystems
    % dimension of the state space, between 2 and 7
    n(i) = 1+unidrnd(nUpperBound-1);
    %number of modes, between 2 and 5
    m(i) = 1+unidrnd(mUpperBound-1);
    A = cell(1,m(i)); %cell A stores the modes of the system
    d = n(i)*(n(i)+1)/2+1;
    
    %number of points in the sampling, between d+1 and 1000
    N(i) = d + unidrnd(NupperBound-1);
    X = cell(1,N(i));
    Y = cell(1,N(i));
    
  
    %generate the m modes
    for j = 1:m(i)
        v=-5+rand(1,n(i))*10;
        D = diag(v); % random eigenvalues of the mode i
        V = orth(randn(n(i)));
        E = V*D*V';
        A{j} = E;
    end
    
    c=jsr_prod_bruteForce(A); %computation of the JSR with brute force method provided by the JSR Toolbox
    upperBound = inf;
    
    while (upperBound == inf)
        [lowerBound, upperBound] = computeRhoBlackbox(A, beta, N(i));
        N(i) = ceil(N(i)*1.5);
    end
    
    fracC(i) = upperBound/c(2);
    if (upperBound >=c(2))
        counterOK=counterOK+1;
        deltaC(counterOK) = (upperBound - c(2));
    else
        deltaC(counterOK) = inf;
    end

end
proba = counterOK/Nsystems;
fprintf('%f\n',proba)
for i = 1:length(m)
    fprintf('The system has %d modes %d states. The number of samples used is: %d\n', m(i), n(i), N(i));
    fprintf('Dfference between the bounds is: %f\n', deltaC(i));
    fprintf('Ratio the bounds is: %f\n', fracC(i));
end