classdef TreeNode < handle
    properties
        parent
        children
        value % node value for mutation
        label % array of cell mutations
        cell_name % if particulat cell has name we can put it here to distinguish later
        freq % cell frequency
        fit % cell fit
        subtree_freq % for non-leaf sum of frequencies of all cells in a subtree
        cell_freq % sum of only dirent children cell frequencies
        node_index
        cell_index
        nodes = [] % just list of all sorted nodes (so we can easily call them)
        cells = [] % same for cells
        time % mutation time
    end
    properties (SetAccess = private, Hidden = true)
        cell_inc = 0 % just utility variable for cloning tree
        pathes = []
        path_index = 1;
        locus_number = 0;
    end
    methods
        %just displays tree
        function dispTree(obj)
            nn = length(obj.nodes);
            nc = length(obj.cells);
            n = nn+nc;
            A(n,n) = 0;
            labels = cell(1,n);
            for i = 1:nn
                if ~isempty(obj.nodes(i).parent)
                    A(obj.nodes(i).parent.node_index, i) = 1;
                end
                labels{i} = [int2str(i) ' m' int2str(obj.nodes(i).value)];
            end
            for i = 1:nc
                A(obj.cells(i).parent.node_index, nn+i) = 1;
                labels{nn+i} = ['c' int2str(i) ' ' num2str(obj.cells(i).freq)];
            end
            bg = biograph(A,labels);
            view(bg);
        end
        %returns copy of the tree with deep copy
        function tree = copy(obj)
            root = TreeNode;
            if ~isempty(obj.nodes)
                n(1,length(obj.nodes))  = TreeNode;
                root.nodes  = n;
            end
            if ~isempty(obj.cells)
                c(1,length(obj.cells)) = TreeNode;
                root.cells = c;
            end
            obj.cell_inc = 1;
            tree = deep_copy(obj, obj, root);
            tree.nodes = root.nodes;
            tree.cells = root.cells;
        end
        
        function copy = deep_copy(obj, node, root)
            copy = TreeNode;
            copy.value = node.value;
            copy.label = node.label;
            copy.freq = node.freq;
            copy.cell_name = node.cell_name;
            copy.node_index = node.node_index;
            copy.cell_index = node.cell_index;
            copy.subtree_freq = node.subtree_freq;
            copy.cell_freq = node.subtree_freq;
            if ~isempty(copy.value)
                root.nodes(copy.node_index) = copy;
            else
                root.cells(obj.cell_inc) = copy;
                obj.cell_inc = obj.cell_inc+1;
            end
            if ~isempty(node.children)
                ch(1,length(node.children)) = TreeNode;
            else
                ch = [];
            end
            copy.children = ch;
            for i = 1:length(node.children)
                child = deep_copy(obj, node.children(i), root);
                copy.children(i) = child;
                child.parent = copy;
            end
        end
        % tells if node has node with given value(mutation) as ancestor
        function has = has_as_ancestor(obj, value)
            has = 0;
            tmp = obj;
            if ~isempty(obj.label)
                tmp = obj.parent;
            end
            while tmp.value ~= 0
                if tmp.value == value
                    has = 1;
                    break;
                end
                tmp = tmp.parent;
            end
        end
        
        % hangs cells to best matching path (in terms of hamming distance)
        % freq is a list of frequencies of each cell
        function infered = hang_cells(obj, cells, freq)
            tmp(size(cells,1)) = TreeNode;
            obj.cells = tmp;
            %path_index is current free index in pathes list
            obj.path_index = 2; % 1 for root path without mutations
            obj.pathes = [];
            obj.locus_number = size(cells,2);
            infered = zeros(size(cells, 1), obj.locus_number);
            all_pathes = 1; % calculate all possible pathes(some are undex X, so we don't count them)
            for i = 1:length(obj.children)
                if ~isempty(obj.children(i).value) && obj.children(i).value ~= 'X'
                    all_pathes = all_pathes + subtree_size(obj, obj.children(i));
                end
            end
            obj.pathes(all_pathes, obj.locus_number+1) = 0;
            obj.pathes(1, obj.locus_number+1) = 1; % n+1 stores index of corresponding node
            for i = 1:length(obj.children)
                if obj.children(i).value ~= 'X'
                    add_path(obj, obj.children(i), obj.pathes(1,1:obj.locus_number));
                end
            end
            % assing cells to nodes with minimum hamming distance
            for i = 1:size(cells,1)
                c = cells(i,:);
                min_k = 1;
                min = sum(xor(c, obj.pathes(1, 1:obj.locus_number)));% hamming distance, assuming that it is bit vectors
                min_node = obj.pathes(1, obj.locus_number + 1);
                for k = 2:size(obj.pathes, 1)
                    d = sum(xor(c, obj.pathes(k, 1:obj.locus_number)));
                    if d < min
                        min = d;
                        min_node = obj.pathes(k, obj.locus_number + 1);
                        min_k = k;
                    end
                end
                infered(i,:) = obj.pathes(min_k, 1:obj.locus_number);
                obj.cells(i).label = cells(i,:);
                obj.cells(i).cell_index = i;
                obj.cells(i).parent = obj.nodes(min_node);
                obj.cells(i).freq = freq(i);
                obj.nodes(min_node).children = [obj.nodes(min_node).children, obj.cells(i)];
            end
        end
        
        function add_path(obj, node, current_path)
            current_path(node.value) = 1;
            obj.pathes(obj.path_index, 1:obj.locus_number) = current_path; % save new path
            obj.pathes(obj.path_index, obj.locus_number+1) = node.node_index; % and index of corresponding node
            obj.path_index = obj.path_index + 1;
            for i = 1:length(node.children)
                if ~isempty(node.children(i).value)
                    add_path(obj, node.children(i), current_path);
                end
            end
        end
        
        function remove_cells(obj)
            obj.cells = [];
            obj.cell_freq = [];
            obj.subtree_freq = [];
            for i = 1:length(obj.children)
                if ~isempty(obj.children(i).label)
                    obj.children(i:end) = [];
                    break;
                else
                   obj.children(i).remove_cells()
                end
            end
            if isempty(obj.children)
                obj.children = [];
            end
        end
        
        % calculates subtree size of a given node(including itself)
        % excluding cells 
        function c = subtree_size(obj, node)
            if ~isempty(node.value)
                c = 1;
                for i = 1:length(node.children)
                    if (~isempty(node.children(i).value))
                        c = c + subtree_size(obj, node.children(i));
                    end
                end
            else
                c = 0;
            end
            
        end
        % hang f node under t if f and t are not directly connected
        % just swap them if f is t parent of t is f parent
        function move_node(obj, f, t)
            from = obj.nodes(f);
            to = obj.nodes(t);
            if to.parent == from
                tmp = to;
                to = from;
                from = tmp;
            end
            if from.parent == to
                tmp = from.value;
                from.value = to.value;
                to.value = tmp;
                return;
            end
            pc = from.parent.children;
            from.parent.children = pc(pc~=from);
            from.parent = to;
            to.children = [to.children from];
        end
        
        function A = get_adjacency_matrix(obj)
            nn = length(obj.nodes);
            nc = length(obj.cells);
            n = nn+nc;
            A(n,n) = 0;
            for i = 1:nn
                if ~isempty(obj.nodes(i).parent)
                    A(obj.nodes(i).parent.value+1, obj.nodes(i).value+1) = 1;
%                     A(obj.nodes(i).value+1,obj.nodes(i).parent.value+1) = 1;
                end
            end
            for i = 1:nc
                A(obj.cells(i).parent.value+1, nn+i) = 1;
%                 A( nn+i, obj.cells(i).parent.value+1) = 1;
            end
        end
    end
end