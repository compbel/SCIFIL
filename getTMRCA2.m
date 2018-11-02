function Tmrca = getTMRCA2(f1,o1,f2,o2,Tmin,Tmax,eps)

% if f1 == f2
%     if eps/(1-eps) == o1/o2
%        Tmrca = Tmin;
%     else
%         Tmrca = Tmax;
%     end
% else
    Tmrca = (log(eps/(1-eps)) + log(o2/o1))/(f2-f1);
    if Tmrca < Tmin
        Tmrca = Tmin;
    end
    if Tmrca > Tmax
        Tmrca = Tmax;
    end
% end