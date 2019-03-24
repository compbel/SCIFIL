function new_stree = transformStree(stree, from, to)
    new_stree = stree;
    f_parent = new_stree{from, 3};
    new_stree{from, 3} = to;
    children = new_stree{f_parent, 1};
    children(children == from) = [];
    new_stree{f_parent, 1} = children;
    new_stree{to, 1} = [new_stree{to, 1} from]; 
end