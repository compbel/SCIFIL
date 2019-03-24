function prob = probMatr(M_infer,M_exper,pFalseNeg,pFalsePos)
n = size(M_infer,1);
 m = size(M_infer,2);

prob = 0;
for i = 1:n
    for j = 1:m
        if M_infer(i,j) == 0
            if M_exper(i,j) == 1
                prob = prob + log(pFalseNeg);
            end
            if M_exper(i,j) == 0
                prob = prob + log(1-pFalseNeg);
            end   
        end
        if M_infer(i,j) == 1
            if M_exper(i,j) == 0
                prob = prob + log(pFalsePos);
            end
            if M_exper(i,j) == 1
                prob = prob + log(1-pFalsePos);
            end   
        end
    end
end