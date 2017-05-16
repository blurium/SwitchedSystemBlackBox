function [A, jsrRaphael] = createRandomSystem(n,m,maxJSR)
%this function takes n and m and creates a random switched linear system
%with maximum JSR upper bound being maxJSR such that non of the individual
%matrices have the same JSR as the JSR of the whole set
jsrRaphael = [0 0];
jsrIndMatrices = [0 0];


while (jsrRaphael(2)-jsrRaphael(1)<0.01) || (jsrRaphael(2) > maxJSR)
for i = 1:m
    sys = drss(n);
    A{i} = 0.9.*sys.a -0.1*randn(n); %maybe this was 0.8
    jsrIndMatrices(i) = max(abs(eig(A{i})));
end
jsrRaphael = jsr_prod_bruteForce(A);
   % jsrIndMatrices
   % jsrRaphael(2)-jsrRaphael(1)
end


end