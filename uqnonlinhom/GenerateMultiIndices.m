function Index = GenerateMultiIndices(Lgerm, order)
N = factorial(Lgerm + order)/(factorial(Lgerm)*factorial(order));
Index = eye(N - 1, Lgerm);
Ind_1 = 1;
Ind_2 = Lgerm;
I = Lgerm;
for p = 2:order
    for L = 1:Lgerm
        TEST = zeros(1, Lgerm);
        TEST(L) = p-1;
        iL = Ind_1;
        while (~isequal(Index(iL,:), TEST))
            iL = iL + 1;
        end
        for k = iL:Ind_2
            I = I + 1;
            Index(I, : ) = Index(L,:) + Index(k,:);
        end
    end
    Ind_1 = Ind_2+1;
    Ind_2 = I;
end
clear I Ind_1 Ind_2 TEST
return